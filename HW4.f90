! Author: Giorgio Torregrosa
! Date: 10/12/2025
! Name: PHZ3151 Homework 4

PROGRAM HW4

IMPLICIT NONE

! Variable declarations 
real :: theta(181), T_s(181), S_theta(181), pressures(3), Ts_old, Ts_new, Ts_current, theta_current, T_a, rho, F_Ts, dF_Ts, a, LHS, RHS, sum, error, d_theta, normalization_factor
real :: p1 = 0.1e5, p2 = 1e5, p3 = 10e5, R_N2 = 296.8, sigma = 5.67e-8, alpha = 0.2, e_a = 0.8, C_p = 1040
real :: C_d = 0.0015, U = 5.0, Newton_error = 0.0001, Balance_error = 1.0, pi = 3.14159265, solution_error = 100, iterations = 0.0, max_iterations = 1000.0, iteration_step = 0.1, acceptable_error = 0.1
integer :: i, j, k, t, p, no_solution_flag = 0

! store pressure values in pressure array
pressures(1) = p1
pressures(2) = p2
pressures(3) = p3 

! initialize step size for later calculations
d_theta = pi / 180.0
normalization_factor = 1 / (2 * pi)

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

! ************************** Begin PART 1A **************************

! open file for writing statistics at each pressure value using netwons method and the rectangle rule
open(unit = 10, file = 'pressure_statistics_rectangle.dat', status = 'replace', action = 'write')

! initialize loop structure for p1
do p = 1, 3

    ! reinitialize variables used in each pressure loop
    iterations = 0
    solution_error = 100
    no_solution_flag = 0

    do while(solution_error > acceptable_error)

        ! Initital guess for T_a, then increment T_a after each iteration
        T_a = 200.0 + iterations

        ! update rho, and calculate turbulent coupling coefficient "a"
        rho = pressures(p) / (R_N2 * T_a)
        a = rho * C_p * C_d * U

        ! inner loop through values of theta, at every value of theta, find the corresponding value of T_s
        do j = 1, 181

            ! initial guess for T_s, reinitialize error
            Ts_old = 250.0
            error = 100.0

            ! use initial guess for T_s to begin newtons method. Iterate through Newtons method to find T_s
            do while(error > Newton_error) 
                F_Ts = (1.0 - alpha) * S_theta(j) + e_a * sigma * T_a**4 - sigma * Ts_old**4 - a * (Ts_old - T_a)
                dF_Ts = -4.0 * sigma * Ts_old**3 - a
                Ts_new = Ts_old - (F_Ts/dF_Ts)
                error = ABS((Ts_new - Ts_old) / Ts_new) * 100.0
                Ts_old = Ts_new
            enddo    

            ! store the found T_s value using Newtons method in the T_s array that is indexed matched with theta array
            T_s(j) = Ts_new
        enddo

        ! after inner loop compare the LHS and RHS of equation 2 with T_s values, if a match is found exit
        ! solve for equation 1 
        LHS = (1.0/4.0) * (1 - alpha) * 866.0

        ! solve for equation 2 using rectangle rule with left endpoints
        sum = 0.0
        do k = 1, 180
            RHS = (((1.0 - e_a) * sigma * T_s(k)**4 + e_a * sigma * T_a**4) * cos(theta(k))) * d_theta
            sum = sum + RHS
        enddo

        ! normalization factor
        sum = sum * normalization_factor

        ! break out of the main loop if an acceptable T_a value has been found
        solution_error = abs(LHS - sum)

        ! increment iterations
        iterations = iterations + iteration_step

        ! check if we've reached max iteration count, if iterations exceed max, break out of loop and flag error for user
        if(iterations > max_iterations) then
            no_solution_flag = 1
            exit
        endif    
    enddo

    ! check whether we converged to a solution, if not, warn the user
    if (no_solution_flag == 0) then
        write(10, *) "For pressure value of:", pressures(p) / 100000
        write(10, *) "Atmospheric temperature found at the above pressure value = ", T_a
        write(10, *) 
        write(10, *) "Theta(degrees)      Surface Temperature"
            do i = 1, 181
                write(10, '(I4, 14X, F12.4)') i - 91, T_s(i)
            enddo
    else
        write(10, *) "NO SOLUTION WAS FOUND, TERMINATING SEARCH..."
    endif
enddo
close(10)
! ************************** END PART 1A ****************************

! ************************** BEGIN PART 1B **************************
!   What are the resultant equatorial (both day- and night-side) and polar surface temperatures for each
!   of the 3 pressure cases? What does this suggest about the effect of pressure on the latitudinal
!   temperature gradient?

! ************************** END PART 1B ****************************

! ************************** BEGIN PART 2A **************************

! open file for writing statistics at each pressure value using netwons method and the midpoint rule
open(unit = 20, file = 'pressure_statistics_midpoint.dat', status = 'replace', action = 'write')

! initialize loop structure for p1
do p = 1, 3

    ! reinitialize variables used in each pressure loop
    iterations = 0
    solution_error = 100
    no_solution_flag = 0

    do while(solution_error > acceptable_error)

        ! Initital guess for T_a, then increment T_a after each iteration
        T_a = 200.0 + iterations

        ! update rho, and calculate turbulent coupling coefficient "a"
        rho = pressures(p) / (R_N2 * T_a)
        a = rho * C_p * C_d * U

        ! inner loop through values of theta, at every value of theta, find the corresponding value of T_s
        do j = 1, 181

            ! initial guess for T_s, reinitialize error
            Ts_old = 250.0
            error = 100.0

            ! use initial guess for T_s to begin newtons method. Iterate through Newtons method to find T_s
            do while(error > Newton_error) 
                F_Ts = (1.0 - alpha) * S_theta(j) + e_a * sigma * T_a**4 - sigma * Ts_old**4 - a * (Ts_old - T_a)
                dF_Ts = -4.0 * sigma * Ts_old**3 - a
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
        LHS = (1.0/4.0) * (1 - alpha) * 866.0

        ! solve for equation 2 using midpoint rule
        sum = 0.0
        do k = 2, 181
            Ts_current = 0.5 * (T_s(k - 1) + T_s(k) )
            theta_current = 0.5 * (theta(k-1) + theta(k))
            RHS = (((1.0 - e_a) * sigma * Ts_current**4 + e_a * sigma * T_a**4) * cos(theta_current)) * d_theta
            sum = sum + RHS
        enddo

        ! normalization factor
        sum = sum * normalization_factor

        ! break out of the main loop if an acceptable T_a value has been found
        solution_error = abs(LHS - sum)

        ! increment iterations
        iterations = iterations + iteration_step

        ! check if we've reached max iteration count, if iterations exceed max, error is flagged for user
        if(iterations > max_iterations) then
            no_solution_flag = 1
            exit
        endif    
    enddo

    ! check whether we converged to a solution, if not, warn the user and end program
    if (no_solution_flag == 0) then
        write(20, *) "For pressure value of:", pressures(p) / 100000
        write(20, *) "Atmospheric temperature found at the above pressure value = ", T_a
        write(20, *) 
        write(20, *) "Theta(degrees)      Surface Temperature"
            do i = 1, 181
                write(20, '(I4, 14X, F12.4)') i - 91, T_s(i)
            enddo
    else
        write(20, *) "NO SOLUTION WAS FOUND, TERMINATING SEARCH..."
    endif
enddo
close(20)
! ************************** END PART 2A ****************************

! ************************** BEGIN PART 2B **************************

! open file for writing statistics at each pressure value using netwons method and the trapezoid rule
open(unit = 30, file = 'pressure_statistics_trapazoid.dat', status = 'replace', action = 'write')

! initialize loop structure for p1
do p = 1, 3

    ! reinitialize variables used in each pressure loop
    iterations = 0
    solution_error = 100
    no_solution_flag = 0

    do while(solution_error > acceptable_error)

        ! Initital guess for T_a, then increment T_a after each iteration
        T_a = 200.0 + iterations

        ! update rho, and calculate turbulent coupling coefficient "a"
        rho = pressures(p) / (R_N2 * T_a)
        a = rho * C_p * C_d * U

        ! inner loop through values of theta, at every value of theta, find the corresponding value of T_s
        do j = 1, 181

            ! initial guess for T_s, reinitialize error
            Ts_old = 250.0
            error = 100.0

            ! use initial guess for T_s to begin newtons method. Iterate through Newtons method to find T_s
            do while(error > Newton_error) 
                F_Ts = (1.0 - alpha) * S_theta(j) + e_a * sigma * T_a**4 - sigma * Ts_old**4 - a * (Ts_old - T_a)
                dF_Ts = -4.0 * sigma * Ts_old**3 - a
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
        LHS = (1.0/4.0) * (1 - alpha) * 866.0

        ! solve for equation 2 using trapezoid rule
        sum = 0.5 * (((1.0 - e_a) * sigma * T_s(1)**4 + e_a * sigma * T_a**4) * cos(theta(1))) * d_theta
        do k = 2, 180
            RHS = (((1.0 - e_a) * sigma * T_s(k)**4 + e_a * sigma * T_a**4) * cos(theta(k))) * d_theta
            sum = sum + RHS
        enddo
        sum = sum + 0.5 * (((1.0 - e_a) * sigma * T_s(181)**4 + e_a * sigma * T_a**4) * cos(theta(181))) * d_theta

        ! normalization factor
        sum = sum * normalization_factor

        ! break out of the main loop if an acceptable T_a value has been found
        solution_error = abs(LHS - sum)

        ! increment iterations
        iterations = iterations + iteration_step

        ! check if we've reached max iteration count, if iterations exceed max, error is flagged for user
        if(iterations > max_iterations) then
            no_solution_flag = 1
            exit
        endif    
    enddo

    ! check whether we converged to a solution, if not, warn the user and end program
    if (no_solution_flag == 0) then
        write(30, *) "For pressure value of:", pressures(p) / 100000
        write(30, *) "Atmospheric temperature found at the above pressure value = ", T_a
        write(30, *) 
        write(30, *) "Theta(degrees)      Surface Temperature"
            do i = 1, 181
                write(30, '(I4, 14X, F12.4)') i - 91, T_s(i)
            enddo
    else
        write(30, *) "NO SOLUTION WAS FOUND, TERMINATING SEARCH..."
    endif
enddo
close(30)
! ************************** END PART 2B ****************************

! ************************** BEGIN PART 2C **************************
! How does the accuracy of the 3 sets of results compare to one another? Why is that?

! ************************** END PART 2C ****************************
END PROGRAM HW4