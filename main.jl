include("inputoutput.jl")

struct Collision
    idx_0
    idx_1
end

function detectCollisionsSimple(grain_data)
    collision_array = Array{Collision,1}()
    for grain_idx = 1:grain_data.num_grains
        for check_grain_idx = grain_idx+1:grain_data.num_grains
            grain_map_idx = (grain_idx-1)*3
            check_grain_map_idx = (check_grain_idx-1)*3
            x_diff = grain_data.q[grain_map_idx+1]-grain_data.q[check_grain_map_idx+1]
            y_diff = grain_data.q[grain_map_idx+2]-grain_data.q[check_grain_map_idx+2]
            dist = sqrt( x_diff*x_diff + y_diff*y_diff )
            if(dist < (grain_data.r[grain_idx]+grain_data.r[check_grain_idx]) )
                #println("Grains $(grain_idx) and $(check_grain_idx) collide")
                #println(" xiff: $(x_diff): ydiff: $(y_diff)")
                #println(" dist: $(dist): r0: $(grain_data.r[grain_idx]) r1: $(grain_data.r[grain_idx])")
                new_collision = Collision(grain_idx,check_grain_idx)
                push!(collision_array,new_collision)
            end
        end
    end
    return collision_array
end

function calculateCollisionForces(grain_data, collision_array)
    collision_forces = 0.0*grain_data.v
    num_of_collisions = length(collision_array)
    return collision_forces

end

function main(args)

    # Read in grain file info
    println(args)
    input_file_name = "testinput.csv"
    if(length(args) > 0)
        if(length(args[1]) == 0)
            input_file_name = args[1]
        else
            println("No input file name specified")
            return
        end
    end
    x,y,theta,r,vx,vy,omega,fixed,rho = readInputCSV(input_file_name)
    q = [x y theta]'[:]
    v = [vx vy omega]'[:]
    # Construct combined mass and moment of inertia vector
    m_linear = pi*r.*r.*rho
    m_angular = 0.5.*m_linear.*r.*r
    m = [m_linear m_linear m_angular]'[:]
    # Create grain data struct
    grain_data = GrainData(q,r,v,fixed,rho,length(x),m)
    println("q: ",grain_data.q)
    println("r: ",grain_data.r)
    println("v: ",grain_data.v)
    println("fixed: ",grain_data.fixed)
    println("rho: ",grain_data.rho)
    println("m: ",grain_data.m)
    #writeOutputVTU("testoutput.vtu",grain_data)

    # Begin time integration
    k = 1000000.0
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
        
        # Forward Euler
        grain_data.v = grain_data.v + (1.0./grain_data.m).*forces*delta_t
        grain_data.q = grain_data.q + grain_data.v*delta_t

        # Handle file output
        ctr += 1
        t += delta_t
        if(ctr%num_of_steps_per_output == 0)
            output_num += 1
            println(ctr," ",t)
            output_file_name = determineOutputFileName("config_",".vtu",max_num_of_outputs,output_num)
            writeOutputVTU(output_file_name,grain_data)
        end
    end
end

main(ARGS)
