include("inputoutput.jl")
include("collisionhandling.jl")

function main(args)
    # Read in grain file info
    println(args)
    #input_file_name = "bouncetest.csv"
    input_file_name = "grainbox.csv"
    #input_file_name = "frictiontest.csv"
    g = [0.0 -9.81]
    #g = [0.0 0.0]
    if(length(args) > 0)
        if(length(args[1]) == 0)
            input_file_name = args[1]
        else
            println("No input file name specified")
            return
        end
    end
    #k_n = 1000000.0
    k_n = 1000000.0
    gamma_n = 75.5575731373*0.5
    k_t = 500000.0
    gamma_t = 37.7787865686*0.5
    mu = 0.4

    x,y,theta,r,vx,vy,omega,fixed,rho = readInputCSV(input_file_name)
    q = [x y theta]'[:]
    v = [vx vy omega]'[:]
    # Construct combined mass and moment of inertia vector
    m_linear = pi*r.*r.*rho
    m_angular = 0.5.*m_linear.*r.*r
    m = [m_linear m_linear m_angular]'[:]
    # Create grain data struct
    indices = 1:length(x)
    indices = reshape(indices,1,length(indices))
    grain_data = GrainData(q,r,v,fixed,rho,length(x),m,k_n,gamma_n,k_t,gamma_t,mu,indices)
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
    delta_t = 0.00001
    output_freq = 1000
    
    num_of_steps_per_sec = Int32(round(1/delta_t))
    max_num_of_outputs = Int32(ceil(output_freq*end_t))
    println("num steps per sec: $(num_of_steps_per_sec)")
    num_of_steps_per_output = Int32(num_of_steps_per_sec/output_freq)
    println("num steps per output: $(num_of_steps_per_output)")
    output_file_name_begin = determineOutputFileName("config_",".vtu",max_num_of_outputs,0)
    writeOutputVTU(output_file_name_begin,grain_data)
    t = 0.0
    output_num = 0
    ctr = 0
    grid_d_mean = 0.0005*2.0
    grid_min = [-4.0*grid_d_mean -4.0*grid_d_mean]
    grid_max = [32.0*grid_d_mean 52.0*grid_d_mean]
    grid_width = grid_d_mean*4.0
    old_collision_array = Array{Collision,1}() 
    while t < end_t
        # Detect collisions
        #@time collisions = detectCollisionsSimple(grain_data,old_collision_array)
        println("=============================================================")
        collisions = detectCollisionsAABB(grain_data,grid_min,grid_max,grid_width,old_collision_array, delta_t)
        println("Collisions: ",length(collisions))
        #println("Collisions:")
        #for ii =1:length(collisions)
        #    println(collisions[ii].idx_0,",",collisions[ii].idx_1)
        #end

        # Calculate collision forces
        println("Forces:")
        @time forces = calculateCollisionForces!(grain_data, collisions,delta_t)

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
        grain_data.v = grain_data.v + forces./grain_data.m*delta_t
        grain_data.q = grain_data.q + grain_data.v*delta_t

        # Handle file output
        ctr += 1
        t += delta_t
        if(ctr%num_of_steps_per_output == 0)
            output_num += 1
            println(ctr," ",t," ",length(collisions))
            output_file_name = determineOutputFileName("config_",".vtu",max_num_of_outputs,output_num)
            writeOutputVTU(output_file_name,grain_data)
            force_output_file_name = determineOutputFileName("force_",".vtu",max_num_of_outputs,output_num)
            writeForcesOutputVTU(force_output_file_name,grain_data,collisions)
        end

        old_collision_array = collisions
    end
end

main(ARGS)
