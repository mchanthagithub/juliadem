include("grainstruct.jl")
using LinearAlgebra

mutable struct Collision
    # always have it so idx_0 < idx_1
    idx_0
    idx_1
    # Collision normal points from idx0 to idx1 
    n # 2 vector, collision normal
    pen_depth # float, penetration depth
    location # 2 vector, collision location
    force # 3 vector, fx, fy, torque
end


# Brute force collision detection
function detectCollisionsSimple(grain_data)
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
                idx_0 = -1
                idx_1 = -1
                if(grain_idx < check_grain_idx)
                    idx_0 = grain_idx
                    idx_1 = check_grain_idx
                else
                    idx_0 = grain_idx
                    idx_1 = check_grain_idx
                end
                map_idx_0 = (idx_0-1)*3
                map_idx_1 = (idx_1-1)*3
                n = [ (grain_data.q[map_idx_1+1] - grain_data.q[map_idx_0+1]) (grain_data.q[map_idx_1+2] - grain_data.q[map_idx_0+2])]
                n = normalize(n)
                pen_depth = (grain_data.r[grain_idx]+grain_data.r[check_grain_idx]) - dist
                location = [0.5*(grain_data.q[map_idx_1+1] + grain_data.q[map_idx_0+1]) 0.5*(grain_data.q[map_idx_1+2] + grain_data.q[map_idx_0+2]) ]
                new_collision = Collision(idx_0,idx_1,n,pen_depth,location,[0.0,0.0,0.0])
                push!(collision_array,new_collision)
            end
        end
    end
    return collision_array
end

function calculateCollisionForces!(grain_data, collision_array)
    collision_forces = 0.0*grain_data.v
    num_of_collisions = length(collision_array)
     
    for collision in collision_array
       idx_0 = collision.idx_0
       idx_1 = collision.idx_1
       map_idx_0 = (idx_0-1)*3
       map_idx_1 = (idx_1-1)*3

       # Normal force is in direction of normal, which goes from idx_0 to idx_1
       normal_force = grain_data.k_n*collision.pen_depth*collision.n

       # Apply collision force to overall force vector
       collision_forces[map_idx_1+1:map_idx_1+2] += normal_force'
       collision_forces[map_idx_0+1:map_idx_0+2] -= normal_force'
       collision.force = [normal_force[1], normal_force[2], 0.0]
    end

    return collision_forces
end

