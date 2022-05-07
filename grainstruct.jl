mutable struct GrainData
    q # 3xn vector, x,y,theta interlaced
    r # n vector, radius
    v # 3xn vector, vx,vy,omega interlaced
    fixed # n vector, whether grain is fixed in place or not
    rho # n vector, density
    num_grains # int, n
    m # 3xn vector, m,m,I interlaced (m is repeated for computational ease)
    k_n # int, normal spring stiffness
end
