! Author: Giorgio Torregrosa
! Date: 10/8/2025
! Name: PHZ3151 Homework 4

PROGRAM HW4

IMPLICIT NONE

! Variable declarations 
real :: theta(181), T_s(181), S_theta(181), d_theta, Ts_old, Ts_new, Ts_current, T_a, rho, F_Ts, dF_Ts, a, Eq_1, Eq_2
real :: p1 = 0.1e5, p2 = 1e5, p3 = 10e5, R_N2 = 296.8, sigma = 5.67e-8, alpha = 0.2, e_a = 0.8, C_p = 1040, C_d = 0.0015, U = 5.0, Newton_error = 0.01, Balance_error = 1.0, pi = 3.14159265
integer :: i, j, iterations = 100


! Begin PART A

! populate theta array
do i = 1, 181
    theta(i) = i - 91
enddo

! populate S_theta array
do i = 1, 181
    S_theta(i) = sin(theta(i) * (pi / 180))*866
enddo        

! initialize loop structure for p1
! outer loop through values of T_a
do i = 0, iterations

    ! increment T_a, update rho, and calculate turbulent coupling coefficient "a"
    T_a = 200 + i
    rho = p1 / (R_N2 * Ta)
    a = rho * C_p * C_d * U

    ! inner loop through values of theta, define values of T_s
    do j = 1, 181

        ! determine if it is night or day
        if(j < 91) then
            F_Ts = e_a * sigma * T_a**4 - sigma * T_s**4 - a*(T_s - T_a)
            dF_Ts = -4*T_s**3 - a
        else 
            F_Ts = (1 - alpha) * 866 * S_theta(j) + e_a * sigma * T_a**4 - sigma * T_s**4 - a*(T_s - T_a)
            dF_Ts = -4*T_s**3 - a
        endif


    enddo
    ! after inner loop compare equations (1) and (2) with T_s values, if a match is found exit

enddo

! repeat the above procedure for the other two different pressure values, 1.0 and 10.0 bars

! END PART A

END PROGRAM HW4