using CSV

mutable struct GrainData
    q
    r
    v
    fixed
    rho
    num_grains
end

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
    fixed = csv_reader.Column8
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
        writeScalerToVTU(file,"Fixed","Int32",grain_data.fixed) 
        writeScalerToVTU(file,"Rho","Float32",grain_data.rho) 
        writeVectorToVTU(file,"Velocity","Float32",grain_data.v,grain_data.num_grains,2,3) 
        writeScalerToVTU(file,"Omega","Float32",grain_data.v[3:3:end]) 
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


x,y,theta,r,vx,vy,omega,fixed,rho = readInputCSV("testinput.csv")
q = [x y theta]'[:]
v = [vx vy omega]'[:]
grain_data = GrainData(q,r,v,fixed,rho,length(x))
println("q: ",grain_data.q)
println("r: ",grain_data.r)
println("v: ",grain_data.v)
println("fixed: ",grain_data.fixed)
println("rho: ",grain_data.rho)

writeOutputVTU("testoutput.vtu",grain_data)
