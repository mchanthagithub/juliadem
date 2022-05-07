include("inputoutput.jl")

struct Collision
    # always have it so idx_0 < idx_1
    idx_0
    idx_1
    # Collision normal points from idx0 to idx1 
    n # 2 vector, collision normal
    pen_depth # float, penetration depth
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
                pen_depth = (grain_data.r[grain_idx]+grain_data.r[check_grain_idx]) - dist
                new_collision = Collision(idx_0,idx_1,n,pen_depth)
                push!(collision_array,new_collision)
            end
        end
    end
    return collision_array
end

function calculateCollisionForces(grain_data, collision_array)
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

    end

    return collision_forces
end

function main(args)

    # Read in grain file info
    println(args)
    #input_file_name = "testinput.csv"
    input_file_name = "bouncetest.csv"
    if(length(args) > 0)
        if(length(args[1]) == 0)
            input_file_name = args[1]
        else
            println("No input file name specified")
            return
        end
    end
    k_n = 1000000.0
    x,y,theta,r,vx,vy,omega,fixed,rho = readInputCSV(input_file_name)
    q = [x y theta]'[:]
    v = [vx vy omega]'[:]
    # Construct combined mass and moment of inertia vector
    m_linear = pi*r.*r.*rho
    m_angular = 0.5.*m_linear.*r.*r
    m = [m_linear m_linear m_angular]'[:]
    # Create grain data struct
    grain_data = GrainData(q,r,v,fixed,rho,length(x),m,k_n)
    println("q: ",grain_data.q)
    println("r: ",grain_data.r)
    println("v: ",grain_data.v)
    println("fixed: ",grain_data.fixed)
    println("rho: ",grain_data.rho)
    println("m: ",grain_data.m)
    println("k: ",grain_data.m)
    #writeOutputVTU("testoutput.vtu",grain_data)

    # Begin time integration
    end_t = 0.1
    delta_t = 0.0001
    output_freq = 1000
    g = [0.0 -9.81]
    num_of_steps_per_sec = Int32(1/delta_t)
    max_num_of_outputs = Int32(ceil(output_freq*end_t))
    println("num steps per sec: $(num_of_steps_per_sec)")
    num_of_steps_per_output = Int32(num_of_steps_per_sec/output_freq)
    println("num steps per output: $(num_of_steps_per_output)")
    output_file_name_begin = determineOutputFileName("config_",".vtu",max_num_of_outputs,0)
    writeOutputVTU(output_file_name_begin,grain_data)
    t = 0.0
    output_num = 0
    ctr = 0
    while t < end_t
        # Detect collisions
        collisions = detectCollisionsSimple(grain_data)
        
        # Calculate collision forces
        forces = calculateCollisionForces(grain_data, collisions)

        # Apply gravity
        f_g_x = grain_data.m[1:3:end].*g[1]
        f_g_y = grain_data.m[2:3:end].*g[2]
        f_g_theta = 0.0*grain_data.m[3:3:end]
        f_g = [f_g_x f_g_y f_g_theta]'[:]
        forces = forces + f_g
       
        # Zero out forces for fixed grains
        for grain = 1:grain_data.num_grains
            if(grain_data.fixed[grain])
                forces[(grain-1)*3+1:(grain-1)*3+3] .= 0
            end
        end

        # Forward Euler
        grain_data.v = grain_data.v + (1.0./grain_data.m).*forces*delta_t
        grain_data.q = grain_data.q + grain_data.v*delta_t

        # Handle file output
        ctr += 1
        t += delta_t
        if(ctr%num_of_steps_per_output == 0)
            output_num += 1
            println(ctr," ",t," ",length(collisions))
            output_file_name = determineOutputFileName("config_",".vtu",max_num_of_outputs,output_num)
            writeOutputVTU(output_file_name,grain_data)
        end
    end
end

main(ARGS)
