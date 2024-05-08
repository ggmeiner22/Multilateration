#=
 * Author: Garrett Gmeiner, ggmeiner2021@my.fit.edu
 * Course: CSE 4250, Fall 2023
 * Project: Project #1, Multilateration
 * Implementation: Julia 1.9.2
=#

using Printf

SPEED_OF_LIGHT = 299792458 * 10^-9  # In Nanometers

satalite_list = []
time_list = []

# Process each line and split the numbers into an array
for line in eachline(stdin)
	numbers_array = parse.(Float64, split(line))
	if length(numbers_array) == 3
        push!(satalite_list, numbers_array)
	else 
	    push!(time_list, numbers_array)
	end
end

# sets the satalites as variables
SATALITE_I = satalite_list[1]
SATALITE_J = satalite_list[2]
SATALITE_K = satalite_list[3]
SATALITE_L = satalite_list[4]

# Preliminary Definitions
x_ji = SATALITE_J[1] - SATALITE_I[1]
x_ki = SATALITE_K[1] - SATALITE_I[1]
x_jk = SATALITE_J[1] - SATALITE_K[1]
x_lk = SATALITE_L[1] - SATALITE_K[1]
y_ki = SATALITE_K[2] - SATALITE_I[2]
y_ji = SATALITE_J[2] - SATALITE_I[2]
y_lk = SATALITE_L[2] - SATALITE_K[2]
y_jk = SATALITE_J[2] - SATALITE_K[2]
z_ji = SATALITE_J[3] - SATALITE_I[3]
z_ki = SATALITE_K[3] - SATALITE_I[3]
z_jk = SATALITE_J[3] - SATALITE_K[3]
z_lk = SATALITE_L[3] - SATALITE_K[3]

# loops through each list of times
for times in time_list
    # More Preliminary Definitions
    r_i = times[1] * SPEED_OF_LIGHT
    r_j = times[2] * SPEED_OF_LIGHT
    r_k = times[3] * SPEED_OF_LIGHT
    r_l = times[4] * SPEED_OF_LIGHT

    # Compute TDOA's
    r_ij = r_i - r_j
    r_ik = r_i - r_k
    r_kl = r_k - r_l
    r_kj = r_k - r_j

    # More Preliminary Definitions
    x_ijy = (r_ij * y_ki) - (r_ik * y_ji)
    x_ikx = (r_ik * x_ji) - (r_ij * x_ki)
    x_ikz = (r_ik * z_ji) - (r_ij * z_ki)
    x_kjy = (r_kj * y_lk) - (r_kl * y_jk)
    x_klx = (r_kl * x_jk) - (r_kj * x_lk)
    x_klz = (r_kl * z_jk) - (r_kj * z_lk)

    # More Preliminary Definitions
    s_i2 = SATALITE_I[1]^2 + SATALITE_I[2]^2 + SATALITE_I[3]^2
    s_j2 = SATALITE_J[1]^2 + SATALITE_J[2]^2 + SATALITE_J[3]^2
    s_k2 = SATALITE_K[1]^2 + SATALITE_K[2]^2 + SATALITE_K[3]^2
    s_l2 = SATALITE_L[1]^2 + SATALITE_L[2]^2 + SATALITE_L[3]^2
 
    # More Preliminary Definitions
    r_ij_2xyz = r_ij^2 + s_i2 - s_j2
    r_ik_2xyz = r_ik^2 + s_i2 - s_k2
    r_kj_2xyz = r_kj^2 + s_k2 - s_j2
    r_kl_2xyz = r_kl^2 + s_k2 - s_l2

    # Saralite Equations
    A = x_ikx / x_ijy
    B = x_ikz / x_ijy
    C = x_klx / x_kjy
    D = x_klz / x_kjy
    
    E = ((r_ik * r_ij_2xyz) - (r_ij * r_ik_2xyz)) / (2 * x_ijy)
    F = ((r_kl * r_kj_2xyz) - (r_kj * r_kl_2xyz)) / (2 * x_kjy)
 
    G = (D - B) / (A - C)
    H = (F - E) / (A - C)

    I = (A * G) + B
    J = (A * H) + E

    K = r_ik_2xyz + (2 * x_ki * H) + (2 * y_ki * J)
    L = 2 * ((x_ki * G) + (y_ki * I) + z_ki)
 
    M = (4 * r_ik^2) * (G^2 + I^2 + 1) - L^2
    N = 8 * r_ik^2 *((G * (SATALITE_I[1] - H)) + (I * (SATALITE_I[2] - J)) + SATALITE_I[3]) + (2 * L * K)
    O = 4 * r_ik^2 * ((SATALITE_I[1] - H)^2 + (SATALITE_I[2] - J)^2 + SATALITE_I[3]^2) - K^2

    Q = N + sign(N) * sqrt(N^2 - (4 * M * O))

    z1 = Q / (2 * M)
    z2 = 2 * O / Q

    Z_pos = 0
    Z_neg = 0
    if sign(z1) == 1
        Z_pos = z1
        Z_neg = z2
    elseif sign(z2) == 1
        Z_pos = z2
        Z_neg = z1
    end

    X_pos = G * Z_pos + H
    Y_pos = I * Z_pos + J

    X_neg = G * Z_neg + H
    Y_neg = I * Z_neg + J

    dis_to_origin_pos = sqrt(X_pos^2 + Y_pos^2 + Z_pos^2)
    dis_to_origin_neg = sqrt(X_neg^2 + Y_neg^2 + Z_neg^2)

    # Using @printf to format and print the text
    @printf("g= %9.2e, h= %9.2e, j= %9.2e, m= %9.2e, o= %9.2e\n", G, H, J, M, O)
    @printf("+) x= %10d, y= %10d, z= %10d; r= %10d\n", X_pos, Y_pos, Z_pos, dis_to_origin_pos)
    @printf("-) x= %10d, y= %10d, z= %10d; r= %10d\n\n", X_neg, Y_neg, Z_neg, dis_to_origin_neg)

end
