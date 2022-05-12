using Random
function main()
# Grain diameter of 1 mm

grain_r_mean = 0.0005
grain_d_mean = grain_r_mean*2.0
#grain_d_mean = 1.0*10^(-3) 
#grain_r_mean = grain_d_mean/2.0

# H = Nd^2/A, where A = Lx*Ly, and Lx and Ly are length and width of bottom plane
# Want H of 40
H = 32
Lx = 32 # grain diameters
N = H*Lx

grain_d_vect = Array{Float64,1}()
grain_r_vect = Array{Float64,1}()
# Generate distribution of diameters
for ii = 1:N
    d = grain_d_mean*2
    # Cutoff so that all grain diameters are +/- 3 sigma of mean (0.85D < D < 1.15D)
    while(d < grain_d_mean*0.85 || d > grain_d_mean*1.15)
      d = grain_d_mean*rand(Float64)
    end
    push!(grain_d_vect,d)
    push!(grain_r_vect,(d/2.0))
end
println(length(grain_r_vect)," ",N)

ctr = 0
spacing = grain_d_mean*1.15*1.1
x = spacing
y = spacing * 2.0
x_coords = Array{Float64,1}()
y_coords = Array{Float64,1}()
while(true)
  ctr = ctr+1
  if(ctr > N)
      d = grain_d_mean*2
      # Cutoff so that all grain diameters are +/- 3 sigma of mean (0.85D < D < 1.15D)
      while(d < grain_d_mean*0.85 || d > grain_d_mean*1.15)
        d = grain_d_mean*rand(Float64)
      end
      push!(grain_d_vect,d)
      push!(grain_r_vect,(d/2.0))
  end
  if(x+grain_r_vect[ctr] > Lx*grain_d_mean-grain_r_mean)
    x = spacing
    y = y + spacing
    if (ctr > N)
      pop!(grain_d_vect)
      pop!(grain_r_vect)
      break
    end
  end
  push!(x_coords,x)
  push!(y_coords,y)
  x = x + spacing
end
#println(x_coords)
#println(y_coords)

ctr = length(x_coords)
N = length(x_coords)
println("Number of actual grains: ",ctr," ",N," ",length(grain_r_vect))

# Material parameters
# Use L3 model from Silbert paper (keep in mind the values from Silber paper are nondimensionalized)
rho = 2600 # Glass beads, kg/m^3
#vol = (4.0/3.0)*math.pi*(grain_r_mean**3)
vol = pi*grain_r_mean*grain_r_mean
m = rho*vol
g = 9.81
#kn = (1.0*10**5)*m*g/grain_d_mean
#kt = (2.0/7.0)*kn
#gamman = 50.0*sqrt(g/grain_d_mean)*m # Don't need to divide by 2 here to get meff because DEM code already multiplies viscous part by 0.50
#gammat = 0.0
mu = 0.40

kn = 1.0e6
kt = kn/2.0
cor = 0.1
gamman = sqrt(m*kn)*(-2.0*log(cor))/sqrt(2.0*(pi*pi + log(cor)*log(cor)))
gammat = gamman/2.0
tauc = pi/sqrt(2.0*kn/m - (gamman/m)*(gamman/m))
println("contact time: ",tauc," dt used: ",string(1.0/250000.0))
# 2.0* is from the fact that in the code gamman is multiplied by 0.5
gamman = gamman * 2.0
gammat = gammat * 2.0
file_name = "grainbox.csv"
open(file_name,"w") do file
    # Create slight velocities to jiggle system a bit
    vel_mean = 0.5*10^(-2)
    println(length(x_coords)," ",length(y_coords)," ",length(grain_r_vect))
    for ctr = 1:N
        write(file,"$(x_coords[ctr])",",","$(y_coords[ctr])",",0.0,","$(grain_r_vect[ctr])",",","$(vel_mean*rand(Float64))",",","$(vel_mean*rand(Float64))",",0.0,0,","$(rho)","\n")  
    end
    ctr = N


    leftCtr = 0
    yLeftCoord = grain_d_mean*H
    xLeftCoord = 0.0 
    while(leftCtr < H)
        write(file,"$(xLeftCoord)",",","$(yLeftCoord)",",0.0,","$(grain_r_mean)",",0.0,0.0,0.0,1,","$(rho)","\n")  
        leftCtr = leftCtr+1
        yLeftCoord = yLeftCoord-grain_d_mean
    end

    botCtr = 0
    yBotCoord = grain_r_mean
    xBotCoord = grain_r_mean
    while(botCtr < Lx)
        write(file,"$(xBotCoord)",",","$(yBotCoord)",",0.0,","$(grain_r_mean)",",0.0,0.0,0.0,1,","$(rho)","\n")  
        botCtr = botCtr+1
        xBotCoord = xBotCoord+grain_d_mean
    end

    rightCtr = 0
    yRightCoord = grain_d_mean*H
    xRightCoord = grain_d_mean*Lx
    while(rightCtr < H)
        write(file,"$(xRightCoord)",",","$(yRightCoord)",",0.0,","$(grain_r_mean)",",0.0,0.0,0.0,1,","$(rho)","\n")  
        rightCtr = rightCtr+1
        yRightCoord = yRightCoord-grain_d_mean
    end

end

end
main()
