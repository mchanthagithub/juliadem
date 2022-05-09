include("grainstruct.jl")
using LinearAlgebra

mutable struct Collision
    # always have it so idx_0 < idx_1
    idx_0 # int, idx of first grain
    idx_1 # int, idx of first grain
    n # 2 vector, collision normal, points from idx0 to idx1 
    pen_depth # float, penetration depth
    location # 2 vector, collision location
    v_rel # 2 vector, velocity of idx0 relative to idx1
    delta_s # 2 vector, tangent spring length
    force # 3 vector, fx, fy, torque
end


# Brute force collision detection
function detectCollisionsSimple(grain_data,old_collision_array)
    collision_array = Array{Collision,1}()
    for grain_idx = 1:grain_data.num_grains
        for check_grain_idx = grain_idx+1:grain_data.num_grains
            grain_map_idx = (grain_idx-1)*3
            check_grain_map_idx = (check_grain_idx-1)*3
            # If both grains are fixed then dont add collision
            if(grain_data.fixed[grain_idx] && grain_data.fixed[check_grain_idx])
                continue
            end

            x_diff = grain_data.q[grain_map_idx+1]-grain_data.q[check_grain_map_idx+1]
            y_diff = grain_data.q[grain_map_idx+2]-grain_data.q[check_grain_map_idx+2]
            dist = sqrt( x_diff*x_diff + y_diff*y_diff )
            if(dist < (grain_data.r[grain_idx]+grain_data.r[check_grain_idx]) )
                #println("Grains $(grain_idx) and $(check_grain_idx) collide")
                #println(" xiff: $(x_diff): ydiff: $(y_diff)")
                #println(" dist: $(dist): r0: $(grain_data.r[grain_idx]) r1: $(grain_data.r[grain_idx])")
                
                # Make sure that idx_0 < idx_1
                idx_0 = -1
                idx_1 = -1
                if(grain_idx < check_grain_idx)
                    idx_0 = grain_idx
                    idx_1 = check_grain_idx
                else
                    idx_1 = grain_idx
                    idx_0 = check_grain_idx
                end
                map_idx_0 = (idx_0-1)*3
                map_idx_1 = (idx_1-1)*3
                n = [ (grain_data.q[map_idx_1+1] - grain_data.q[map_idx_0+1]) (grain_data.q[map_idx_1+2] - grain_data.q[map_idx_0+2])]
                n = normalize(n)
                pen_depth = (grain_data.r[grain_idx]+grain_data.r[check_grain_idx]) - dist
                location = [0.5*(grain_data.q[map_idx_1+1] + grain_data.q[map_idx_0+1]) 0.5*(grain_data.q[map_idx_1+2] + grain_data.q[map_idx_0+2]) ]
                v_0_angular_loc = location - grain_data.q[map_idx_0+1:map_idx_0+2]
                v_0_angular = grain_data.v[map_idx_0+3].*[-v_0_angular_loc[2] v_0_angular_loc[1]]
                v_0_contact = grain_data.v[map_idx_0+1:map_idx_0+2] + v_0_angular 
                v_1_angular_loc = location - grain_data.q[map_idx_1+1:map_idx_1+2]
                v_1_angular = grain_data.v[map_idx_1+3].*[-v_1_angular_loc[2] v_1_angular_loc[1]]
                v_1_contact = grain_data.v[map_idx_1+1:map_idx_1+2] + v_1_angular
                #v_rel = [ (grain_data.v[map_idx_0+1] - grain_data.v[map_idx_1+1]) (grain_data.v[map_idx_0+2] - grain_data.v[map_idx_1+2])]
                v_rel = v_0_contact - v_1_contact 
                
                delta_s = [0.0 0.0]
                for ii = 1:length(old_collision_array)
                    if ( ((old_collision_array[ii].idx_0 == idx_0) && (old_collision_array[ii].idx_1 == idx_1)) || 
                        ((old_collision_array[ii].idx_0 == idx_1) && (old_collision_array[ii].idx_1 == idx_0)) )
                        delta_s = old_collision_array[ii].delta_s + dt*v_rel
                        delta_s -= dot(n,delta_s)*n
                        break
                    end
                end
                
                new_collision = Collision(idx_0,idx_1,n,pen_depth,location,v_rel,delta_s,[0.0,0.0,0.0])
                push!(collision_array,new_collision)
            end
        end
    end
    return collision_array
end

function generateAABBs(grain_data)
    aabb_mins = Matrix{Float64}[]
    aabb_maxes = Matrix{Float64}[]
    aabb_grain_indices = grain_data.indices

    for grain_num = 1:grain_data.num_grains
        grain_x = grain_data.q[(grain_num-1)*3+1]
        grain_y = grain_data.q[(grain_num-1)*3+2]
        grain_r = grain_data.r[grain_num]

        aabb_grain_min = [ (grain_x-grain_r) (grain_y-grain_r)]
        aabb_grain_max = [ (grain_x+grain_r) (grain_y+grain_r)]
        push!(aabb_mins,aabb_grain_min)
        push!(aabb_maxes,aabb_grain_max)
    end
    return aabb_mins,aabb_maxes,aabb_grain_indices
end

function gridAABBProjection!(aabb_mins,aabb_maxes,aabb_grain_indices,grid)
    num_aabbs = length(aabb_mins)
    for aabb_num = 1:num_aabbs
        #println("min: ",aabb_num,": ",(aabb_mins[aabb_num] - grid.min)./grid.width)
        #println("max: ",aabb_num,": ",(aabb_maxes[aabb_num] - grid.min)./grid.width)
        grid_idx_collision_min = floor.(Int32,(aabb_mins[aabb_num] - grid.min)./grid.width).+1
        grid_idx_collision_max = ceil.(Int32,(aabb_maxes[aabb_num] - grid.min)./grid.width).+1
        for kk = grid_idx_collision_min[2]:grid_idx_collision_max[2]
            for ii = grid_idx_collision_min[1]:grid_idx_collision_max[1]
                grid_flat_idx = ii + (kk-1)*grid.num_elems[1]
                #grid.aabbs_in_elem[grid_flat_idx]
                #aabb_grain_indices[aabb_num]
                push!(grid.aabbs_in_elem[grid_flat_idx],aabb_grain_indices[aabb_num])
            end
        end
    end
end

function possibleAABBCollisions(aabb_mins,aabb_maxes,aabb_grain_indices,grid)
    num_grid_elems = length(grid.aabbs_in_elem)
    possible_collisions = Set{Array{Int32, 1}}()


    collision_found = false
    for elem = 1:num_grid_elems
        num_aabbs_in_elem = length(grid.aabbs_in_elem[elem])
        if(num_aabbs_in_elem == 0)
            continue
        end

        for aabb_num = 1:num_aabbs_in_elem
            aabb_actual_grain_idx = grid.aabbs_in_elem[elem][aabb_num]
            for aabb_check_num = aabb_num+1:num_aabbs_in_elem
                aabb_check_actual_grain_idx = grid.aabbs_in_elem[elem][aabb_check_num]
                if(aabb_mins[aabb_actual_grain_idx][1] < aabb_maxes[aabb_check_actual_grain_idx][1] &&
                   aabb_mins[aabb_actual_grain_idx][2] < aabb_maxes[aabb_check_actual_grain_idx][2] &&
                   aabb_maxes[aabb_actual_grain_idx][1] > aabb_mins[aabb_check_actual_grain_idx][1] &&
                   aabb_maxes[aabb_actual_grain_idx][2] > aabb_mins[aabb_check_actual_grain_idx][2])
                   
                    # Make sure that idx_0 < idx_1
                    idx_0 = -1
                    idx_1 = -1
                    if(aabb_actual_grain_idx < aabb_check_actual_grain_idx)
                        idx_0 = aabb_actual_grain_idx
                        idx_1 = aabb_check_actual_grain_idx
                    else
                        idx_1 = aabb_actual_grain_idx
                        idx_0 = aabb_check_actual_grain_idx
                    end
                    
                    push!(possible_collisions,[idx_0,idx_1])
                end
            end
        end
    end

    #if(collision_found)
    #    for elem = 1:num_grid_elems
    #        print("Elem: ",elem,": ")
    #        for ii = 1:length(grid.aabbs_in_elem[elem])
    #            print(grid.aabbs_in_elem[elem][ii]," ")
    #        end
    #        print("\n")
    #    end
    #end
    #println(possible_collisions)

    return possible_collisions

end

function finalStateCollisionDetection(grain_data,possible_collisions,old_collision_array, dt) 
    collision_array = Array{Collision,1}()
    num_of_possible_collisions = length(possible_collisions)
    
    for collision_check in possible_collisions 
        #grain_idx = possible_collisions[possible_collision_num][1]
        #check_grain_idx = possible_collisions[possible_collision_num][2]
        grain_idx = collision_check[1]
        check_grain_idx = collision_check[2]

        # If both grains are fixed then dont add collision
        if(grain_data.fixed[grain_idx] && grain_data.fixed[check_grain_idx])
            continue
        end
        grain_map_idx = (grain_idx-1)*3
        check_grain_map_idx = (check_grain_idx-1)*3
        x_diff = grain_data.q[grain_map_idx+1]-grain_data.q[check_grain_map_idx+1]
        y_diff = grain_data.q[grain_map_idx+2]-grain_data.q[check_grain_map_idx+2]
        dist = sqrt( x_diff*x_diff + y_diff*y_diff )
        if(dist < (grain_data.r[grain_idx]+grain_data.r[check_grain_idx]) )
            #println("Grains $(grain_idx) and $(check_grain_idx) collide")
            #println(" xiff: $(x_diff): ydiff: $(y_diff)")
            #println(" dist: $(dist): r0: $(grain_data.r[grain_idx]) r1: $(grain_data.r[grain_idx])")
            
            # Make sure that idx_0 < idx_1
            idx_0 = -1
            idx_1 = -1
            if(grain_idx < check_grain_idx)
                idx_0 = grain_idx
                idx_1 = check_grain_idx
            else
                idx_1 = grain_idx
                idx_0 = check_grain_idx
            end
            map_idx_0 = (idx_0-1)*3
            map_idx_1 = (idx_1-1)*3
            n = [ (grain_data.q[map_idx_1+1] - grain_data.q[map_idx_0+1]) (grain_data.q[map_idx_1+2] - grain_data.q[map_idx_0+2])]
            n = normalize(n)
            pen_depth = (grain_data.r[grain_idx]+grain_data.r[check_grain_idx]) - dist
            location = [0.5*(grain_data.q[map_idx_1+1] + grain_data.q[map_idx_0+1]) 0.5*(grain_data.q[map_idx_1+2] + grain_data.q[map_idx_0+2]) ]
            v_0_angular_loc = location - [grain_data.q[map_idx_0+1] grain_data.q[map_idx_0+2]] 
            v_0_angular = grain_data.v[map_idx_0+3].*[-v_0_angular_loc[2] v_0_angular_loc[1]]
            v_0_contact = [grain_data.v[map_idx_0+1] grain_data.v[map_idx_0+2]]+ v_0_angular 
            v_1_angular_loc = location - [grain_data.q[map_idx_1+1] grain_data.q[map_idx_1+2]] 
            v_1_angular = grain_data.v[map_idx_1+3].*[-v_1_angular_loc[2] v_1_angular_loc[1]]
            v_1_contact = [grain_data.v[map_idx_1+1] grain_data.v[map_idx_1+2]] + v_1_angular
            #v_rel = [ (grain_data.v[map_idx_0+1] - grain_data.v[map_idx_1+1]) (grain_data.v[map_idx_0+2] - grain_data.v[map_idx_1+2])]
            v_rel = v_0_contact - v_1_contact 
            v_rel_old = [ (grain_data.v[map_idx_0+1] - grain_data.v[map_idx_1+1]) (grain_data.v[map_idx_0+2] - grain_data.v[map_idx_1+2])]
            #println(v_rel)
            #println(v_rel_old)
            
            delta_s = [0.0 0.0]
            found_old_spring = false
            for ii = 1:length(old_collision_array)
                if ( ((old_collision_array[ii].idx_0 == idx_0) && (old_collision_array[ii].idx_1 == idx_1)) || 
                    ((old_collision_array[ii].idx_0 == idx_1) && (old_collision_array[ii].idx_1 == idx_0)) )
                    delta_s = old_collision_array[ii].delta_s + dt*v_rel
                    delta_s -= dot(n,delta_s)*n
                    found_old_spring = true
                    break
                end
            end

            if (!found_old_spring)
                delta_s = dt*v_rel
                delta_s -= dot(n,delta_s)*n
            end
            
            new_collision = Collision(idx_0,idx_1,n,pen_depth,location,v_rel,delta_s,[0.0,0.0,0.0])
            #println(new_collision)
            push!(collision_array,new_collision)
        end
    end
    
    return collision_array


end

mutable struct CollisionGrid
    min
    max
    width
    num_elems
    aabbs_in_elem
end

function initAABBVector(num_grid_elems)
    aabbs_in_elem = Vector{Vector{Int32}}(undef,0)
    for ii = 1:num_grid_elems[1]*num_grid_elems[2]
        temp = Vector{Int32}(undef,0)
        push!(aabbs_in_elem,temp)
    end
    return aabbs_in_elem
end

# Collision detection using a grid and axis-aligned bounding boxes
function detectCollisionsAABB(grain_data,grid_min,grid_max,grid_width, old_collision_array, dt)
    println("generate aabb: ")
    @time aabb_mins,aabb_maxes,aabb_grain_indices = generateAABBs(grain_data)
    #aabbs_in_elem = Vector{Vector{Int32}}(undef,0)
    num_grid_elems = ((grid_max.+0.00001).-grid_min)./grid_width
    num_grid_elems = floor.(Int32,num_grid_elems)
    println("init aabb_in_elem: ")
    @time aabbs_in_elem = initAABBVector(num_grid_elems) 
    println("init collision: ")
    @time collision_grid = CollisionGrid(grid_min,grid_max,grid_width,num_grid_elems,aabbs_in_elem)
    println("project aabb: ")
    @time gridAABBProjection!(aabb_mins,aabb_maxes,aabb_grain_indices,collision_grid)
    println("possible collision: ")
    @time possible_collisions = possibleAABBCollisions(aabb_mins,aabb_maxes,aabb_grain_indices,collision_grid)
    println("final collision: ")
    @time collision_array = finalStateCollisionDetection(grain_data,possible_collisions, old_collision_array, dt)
    return collision_array
end

function secondRootOfQuadratic(a,b,c,dscr_sqrt)
    if(b > 0)
        root = 2*c / (-b - dscr_sqrt)
    else
        root = (-b + dscr_sqrt) / (2.0 * a) 
    end
    return root
end


function calculateCollisionForces!(grain_data, collision_array, dt)
    collision_forces = 0.0*grain_data.v
    num_of_collisions = length(collision_array)
     
    for collision in collision_array
       idx_0 = collision.idx_0
       idx_1 = collision.idx_1
       map_idx_0 = (idx_0-1)*3
       map_idx_1 = (idx_1-1)*3
       #println("Coll: ",collision) 
       #println("idx0: ",idx_0)
       #println("idx1: ",idx_0)
       # Normal force is in direction of normal, which goes from idx_0 to idx_1
       #println("Spring")
       #println(grain_data.k_n*collision.pen_depth*collision.n)
       #println("Damper")
       #println(grain_data.gamma_n*collision.v_rel)

       v_rel_tangent = collision.v_rel - dot(collision.n,collision.v_rel)*collision.n
       v_rel_n = collision.v_rel - v_rel_tangent
       #println("vreltan: ",v_rel_tangent)
       normal_force = grain_data.k_n*collision.pen_depth*collision.n .+ grain_data.gamma_n*v_rel_n
       friction_force = -grain_data.k_t*collision.delta_s .- grain_data.gamma_n*v_rel_tangent
       mu_fn = grain_data.mu * norm(normal_force)
       ft = norm(friction_force)
       #println("friction force: ",friction_force)
       if( 0.5 * grain_data.gamma_t * norm(v_rel_tangent) > mu_fn )
           collision.delta_s = [0.0 0.0]
           friction_force = - grain_data.gamma_t*v_rel_tangent
           friction_force = normalize(friction_force)
           friction_force *= mu_fn
           #println("friction force a: ",friction_force)
       elseif( ft > mu_fn )
           a = grain_data.k_t * grain_data.k_t * norm(collision.delta_s)*norm(collision.delta_s)
           #println("a: ",a)
           b = grain_data.k_t * grain_data.gamma_t * dot(collision.delta_s, v_rel_tangent )
           #println("b: ",b)
           c = 0.25 * grain_data.gamma_t * grain_data.gamma_t * dot(v_rel_tangent,v_rel_tangent) - grain_data.mu * grain_data.mu * norm(normal_force)*norm(normal_force)
           #println("c: ",c)
           dscr =  b * b - 4 * a * c 
           #println("dscr: ",dscr)
           dscr_sqrt = sqrt(dscr)
           #println("dscr_sqrt: ",dscr_sqrt)
           root1 = secondRootOfQuadratic( a, b, c, dscr_sqrt )
           #println("root1: ",root1)
           new_delta_s =  root1 * collision.delta_s
           #println("news: ",new_delta_s)
           new_friction_force = -grain_data.k_t * new_delta_s - grain_data.gamma_t * v_rel_tangent 
           #println("newf: ",new_friction_force)
           friction_force = new_friction_force
           collision.delta_s = new_delta_s
           #println("friction force b: ",friction_force)
       end
       friction_force = -friction_force
       
       #println("friction force: ",friction_force)
       total_force = normal_force + friction_force
       #println("total force: ",total_force)
       # Apply collision force to overall force vector
       collision_forces[map_idx_1+1:map_idx_1+2] += total_force'
       collision_forces[map_idx_0+1:map_idx_0+2] -= total_force'
      
       lever_1 = collision.location - [grain_data.q[map_idx_1+1] grain_data.q[map_idx_1+2]]
       lever_0 = collision.location - [grain_data.q[map_idx_0+1] grain_data.q[map_idx_0+2]]
       #println("lever1: ",lever_1)
       #println("lever0: ",lever_0)
       torque_1 = lever_1[1]*total_force[2] - lever_1[2]*total_force[1] 
       torque_0 = lever_0[1]*total_force[2] - lever_0[2]*total_force[1]
       #println("torque_1: ",torque_1)
       #println("torque_0: ",torque_0)
       collision_forces[map_idx_1+3] += torque_1
       collision_forces[map_idx_0+3] -= torque_0

       # Save out the collision force into the collision data structure
       collision.force = [normal_force[1], normal_force[2], 0.0]
    end

    return collision_forces
end

