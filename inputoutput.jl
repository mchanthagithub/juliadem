using CSV
include("grainstruct.jl")

function readInputCSV(file_name)
    csv_reader = CSV.File(file_name,header=0)
    #for row in csv_reader
    #    println(row)
    #end
    #println(csv_reader.names)
    x = csv_reader.Column1
    y = csv_reader.Column2
    theta = csv_reader.Column3
    r = csv_reader.Column4
    vx = csv_reader.Column5
    vy = csv_reader.Column6
    omega = csv_reader.Column7
    fixed = !=(0).(csv_reader.Column8) # Need to convert the 1 and 0 to bools
    rho = csv_reader.Column9

    return x,y,theta,r,vx,vy,omega,fixed,rho
end

function writeScalerToVTU(file,field_name,field_type,field_data)
    write(file,"    <DataArray type=\"$(field_type)\" Name=\"$(field_name)\" format=\"ascii\">\n")
    ctr = 0
    spaceFlag = false
    width = 10
    for kk in 1:length(field_data)
        spaceFlag = false
        write(file,"$(field_data[kk]) ")
        ctr = ctr + 1
        if(ctr%width == 0)
            write(file,"\n")
            spaceFlag = true    
        end 
    end
    if(!spaceFlag)
        write(file,"\n")
    end
    write(file,"    </DataArray>\n")
end

function writeVectorToVTU(file,field_name,field_type,field_data,num_data_pts,num_to_extract,data_offset)
    write(file,"    <DataArray type=\"$(field_type)\" Name=\"$(field_name)\" NumberOfComponents=\"3\" format=\"ascii\">\n")
    ctr = 0
    spaceFlag = false
    width = 10
    for kk in 1:num_data_pts
        spaceFlag = false
        write(file,"$(field_data[(kk-1)*data_offset+1]) ")
        for ii in 2:num_to_extract
            write(file,"$(field_data[(kk-1)*data_offset+ii]) ")
        end
        if(num_to_extract < 3)
            write(file,"0.0 ")
        end
        ctr = ctr + 1
        if(ctr%width == 0)
            write(file,"\n")
            spaceFlag = true    
        end 
    end
    if(!spaceFlag)
        write(file,"\n")
    end
    write(file,"    </DataArray>\n")
end

function write2DVectorToVTU(file,field_name,field_type,field_data,num_data_pts)
    write(file,"    <DataArray type=\"$(field_type)\" Name=\"$(field_name)\" NumberOfComponents=\"2\" format=\"ascii\">\n")
    ctr = 0
    spaceFlag = false
    width = 10
    for kk in 1:num_data_pts
        spaceFlag = false
        write(file,"$(field_data[(kk-1)*2+1]) ")
        write(file,"$(field_data[(kk-1)*2+2]) ")
        ctr = ctr + 1
        if(ctr%width == 0)
            write(file,"\n")
            spaceFlag = true    
        end 
    end
    if(!spaceFlag)
        write(file,"\n")
    end
    write(file,"    </DataArray>\n")
end

function writeOutputVTU(file_name,grain_data)
    open(file_name,"w") do file
        write(file,"<?xml version=\"1.0\"?>\n")
        write(file,"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n")
        write(file," <UnstructuredGrid>\n")
        write(file,"  <Piece NumberOfPoints=\"$(grain_data.num_grains)\" NumberOfCells=\"0\">\n")
        write(file,"   <Points>\n")
        write(file,"    <DataArray name=\"Position\" type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n")
        for kk in 1:grain_data.num_grains
             write(file,"$(grain_data.q[(kk-1)*3+1]) $(grain_data.q[(kk-1)*3+2]) 0.0\n")
        end
        write(file,"    </DataArray>\n")
        write(file,"   </Points>\n")
        write(file,"   <PointData Scalars=\"scalars\">\n")
        writeVectorToVTU(file,"Position","Float32",grain_data.q,grain_data.num_grains,2,3) 
        writeScalerToVTU(file,"Theta","Float32",grain_data.q[3:3:end]) 
        writeScalerToVTU(file,"Radius","Float32",grain_data.r) 
        writeScalerToVTU(file,"Fixed","Int32",Int32[x ? 1 : 0 for x in grain_data.fixed]) # Need to convert bool to ints
        writeScalerToVTU(file,"Rho","Float32",grain_data.rho) 
        writeVectorToVTU(file,"Velocity","Float32",grain_data.v,grain_data.num_grains,2,3) 
        writeScalerToVTU(file,"Omega","Float32",grain_data.v[3:3:end]) 
        writeVectorToVTU(file,"Mass","Float32",grain_data.m,grain_data.num_grains,2,3) 
        writeScalerToVTU(file,"MomentOfInertia","Float32",grain_data.m[3:3:end]) 
        writeScalerToVTU(file,"Indices","Int32",grain_data.indices) 
        write(file,"   </PointData>\n")
        write(file,"   <Cells>\n")
        write(file,"    <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n")
        write(file,"    </DataArray>\n")
        write(file,"    <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n")
        write(file,"    </DataArray>\n")
        write(file,"    <DataArray type=\"Int32\" Name=\"types\" format=\"ascii\">\n")
        write(file,"    </DataArray>\n")
        write(file,"   </Cells>\n")
        write(file,"  </Piece>\n")
        write(file," </UnstructuredGrid>\n")
        write(file,"</VTKFile>\n")
    end
end

function writeForcesOutputVTU(file_name,grain_data,collision_data_array)
    num_of_collisions = length(collision_data_array)
    idx_0 = Array{Int32,1}()
    idx_1 = Array{Int32,1}()
    n = Array{Float64,1}()
    pen_depth = Array{Float64,1}()
    location = Array{Float64,1}()
    force = Array{Float64,1}()

    if(num_of_collisions == 0)
        push!(idx_0,-1)
        push!(idx_1,-1)
        push!(n,0.0)
        push!(n,0.0)
        push!(pen_depth,0.0)
        push!(location,0.0)
        push!(location,0.0)
        push!(force,0.0)
        push!(force,0.0)
        push!(force,0.0)
    end

    for ii = 1:num_of_collisions
        push!(idx_0,collision_data_array[ii].idx_0)
        push!(idx_1,collision_data_array[ii].idx_1)
        push!(n,collision_data_array[ii].n[1])
        push!(n,collision_data_array[ii].n[2])
        push!(pen_depth,collision_data_array[ii].pen_depth)
        push!(location,collision_data_array[ii].location[1])
        push!(location,collision_data_array[ii].location[2])
        push!(force,collision_data_array[ii].force[1])
        push!(force,collision_data_array[ii].force[2])
        push!(force,0.0)
    end

    if(num_of_collisions == 0)
        num_of_collisions = 1
    end

    open(file_name,"w") do file
        write(file,"<?xml version=\"1.0\"?>\n")
        write(file,"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n")
        write(file," <UnstructuredGrid>\n")
        write(file,"  <Piece NumberOfPoints=\"$(num_of_collisions)\" NumberOfCells=\"0\">\n")
        write(file,"   <Points>\n")
        write(file,"    <DataArray name=\"Position\" type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n")
        for kk in 1:num_of_collisions
             write(file,"$(location[(kk-1)*2+1]) $(location[(kk-1)*2+2]) 0.0\n")
        end
        write(file,"    </DataArray>\n")
        write(file,"   </Points>\n")
        write(file,"   <PointData Scalars=\"scalars\">\n")
        write2DVectorToVTU(file,"Position","Float32",location,num_of_collisions) 
        writeScalerToVTU(file,"Idx0","Int32",idx_0) 
        writeScalerToVTU(file,"Idx1","Int32",idx_1) 
        write2DVectorToVTU(file,"n","Float32",n,num_of_collisions) 
        writeScalerToVTU(file,"PenDepth","Float32",pen_depth) 
        writeVectorToVTU(file,"Force","Float32",force,num_of_collisions,3,3) 
        write(file,"   </PointData>\n")
        write(file,"   <Cells>\n")
        write(file,"    <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n")
        write(file,"    </DataArray>\n")
        write(file,"    <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n")
        write(file,"    </DataArray>\n")
        write(file,"    <DataArray type=\"Int32\" Name=\"types\" format=\"ascii\">\n")
        write(file,"    </DataArray>\n")
        write(file,"   </Cells>\n")
        write(file,"  </Piece>\n")
        write(file," </UnstructuredGrid>\n")
        write(file,"</VTKFile>\n")
    end
end



function determineOutputFileName(prefix,postfix,max_num_of_outputs,output_num)
    max_num_of_digits = Int32(floor(log10(max_num_of_outputs)+1))
    if(output_num == 0)
        output_num_of_digits = 1
    else
        output_num_of_digits = Int32(floor(log10(output_num)+1))
    end
    num_of_leading_zeros = max_num_of_digits - output_num_of_digits
    file_name = prefix*(("0")^num_of_leading_zeros)*string(output_num)*postfix
    #println("max num digits: $(max_num_of_digits) num: $(max_num_of_outputs)")
    #println("output num digits: $(output_num_of_digits) num: $(output_num)")
    println(file_name)
    return file_name
end
