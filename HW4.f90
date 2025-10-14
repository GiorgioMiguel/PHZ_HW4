! Author: Giorgio Torregrosa
! Date: 10/12/2025
! Name: PHZ3151 Homework 4

PROGRAM HW4

IMPLICIT NONE

! Variable declarations 
real :: theta(181), T_s(181), S_theta(181), d_theta, Ts_old, Ts_new, Ts_current, T_a, rho, F_Ts, dF_Ts, a, LHS, RHS, sum, error
real :: p1 = 0.1e5, p2 = 1e5, p3 = 10e5, R_N2 = 296.8, sigma = 5.67e-8, alpha = 0.2, e_a = 0.8, C_p = 1040, C_d = 0.0015, U = 5.0, Newton_error = 0.0001, Balance_error = 1.0, pi = 3.14159265, solution_error = 100
integer :: i, j, k, t, iterations = 0, max_iterations = 10000, no_solution_flag = 0

! Begin PART A

! populate theta array with radian values
do i = 1, 181
    theta(i) = (i - 91) * (pi / 180.0) 
enddo

! populate S_theta array
do i = 1, 181
    if(i - 91 <= 0) then
        S_theta(i) = 0
    else
        S_theta(i) = sin(theta(i)) * 866.0
    endif    
enddo        

! initialize loop structure for p1
do while(solution_error > 1.0)

    ! Initital guess for T_a, then increment T_a after each iteration
    T_a = 200.0 + iterations

    ! update rho, and calculate turbulent coupling coefficient "a"
    rho = p1 / (R_N2 * T_a)
    a = rho * C_p * C_d * U

    ! inner loop through values of theta, at every value of theta, find the corresponding value of T_s
    do j = 1, 181

        ! initial guess for T_s, reinitialize error
        Ts_old = 250.0
        error = 100.0

        ! use initial guess for T_s to begin newtons method. Iterate through Newtons method to find T_s
        do while(error > Newton_error) 
            F_Ts = (1.0 - alpha) * S_theta(j) + e_a * sigma * T_a**4 - sigma * Ts_old**4 - a*(Ts_old - T_a)
            dF_Ts = -4.0*sigma*Ts_old**3 - a
            Ts_new = Ts_old - (F_Ts/dF_Ts)
            error = ABS((Ts_new - Ts_old) / Ts_new) * 100.0
            Ts_old = Ts_new
        enddo    

        ! store the found T_s value using Newtons method in the T_s array that is indexed matched with theta array
        T_s(j) = Ts_new
    enddo

    ! after inner loop compare equations (1) and (2) with T_s values, if a match is found exit
    ! initialize d_theta step size
    d_theta = pi / 180.0
    
    ! solve for equation 1 
    LHS = (1.0/4.0)*(1 - alpha) * 866.0

    ! solve for equation 2 using rectangle rule
    sum = 0.0
    do k = 1, 180
        RHS = (((1.0 - e_a)*sigma*T_s(k)**4 + e_a * sigma * T_a**4) * cos(theta(k))) * d_theta
        sum = sum + RHS
    enddo

    ! normalization factor
    sum = sum * (1/(2*pi))

    ! break out of the main loop if an acceptable T_a value has been found
    solution_error = abs(LHS - sum)

    ! increment iterations
    iterations = iterations + 1

    ! check if we've reached max iteration count, if iterations exceed max, error is flagged for user
    if(iterations > max_iterations) then
        no_solution_flag = 1
        exit
    endif    
enddo

! check whether we converged to a solution, if not, warn the user and end program
if (no_solution_flag == 0) then
    open(unit = 10, file = '0.1_pressure_statistics.dat', status = 'replace', action = 'write')
    write(10, *) "Atmospheric temperature found at pressure value of 0.1 bars = ", T_a
    write(10, *) 
    write(10, *) "Theta(degrees)      Surface Temperature"
        do i = 1, 181
            write(10, '(I4, 14X, F12.4)') i-91, T_s(i)
        enddo
    close(10)
else
    print *, "NO SOLUTION WAS FOUND, TERMINATING PROGRAM..."
endif

! repeat the above procedure for the other two different pressure values, 1.0 and 10.0 bars
! (NOT IMPLEMNTED YET)

! END PART A

END PROGRAM HW4