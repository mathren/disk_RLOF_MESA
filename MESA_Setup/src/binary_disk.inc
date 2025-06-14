! TODO: this uses Kramers opacity law for simplicity, switch to
! tabulated opacity for self-consistency

! ----------------------------------------------------------------------------------------------
subroutine eval_L2_mass_loss_fraction(donor_mass, accretor_mass, mass_transfer_rate, orbital_separation, &
     disk_alpha, disk_mu, &
     fL2, ierr)

  ! Calculate the (outer) L2 mass-loss fraction
  ! according to Lu et al. (2023, MNRAS 519, 1409) "On rapid binary mass transfer -I. Physical model"
  ! https://ui.adsabs.harvard.edu/abs/2023MNRAS.519.1409L/abstract

  ! constants
  use const_def, only: dp, pi, pi2, one_third, two_thirds, standard_cgrav, Msun, Rsun, secyer, crad, boltzm, clight, mp
  real(dp), parameter :: one_fourth = 1.0_dp / 4.0_dp
  real(dp), parameter :: one_eigth = 1.0_dp / 8.0_dp
  real(dp), parameter :: three_eigths = 3.0_dp / 8.0_dp
  real(dp), parameter :: one_ninth = 1.0_dp / 9.0_dp
  real(dp), parameter :: four_twentyseventh = 4.0_dp / 27.0_dp

  real(dp), intent(in) :: donor_mass         ! [M_sun]
  real(dp), intent(in) :: accretor_mass      ! [M_sun]
  real(dp), intent(in) :: mass_transfer_rate ! [M_sun/yr]
  real(dp), intent(in) :: orbital_separation ! [R_sun]
  real(dp), intent(in) :: disk_alpha         ! disk alpha viscosity parameter (dimensionless)
  real(dp), intent(in) :: disk_mu            ! disk mean molecular weight (dimensionless)
  real(dp), intent(out) :: fL2               ! L2 mass-loss fraction (dimensionless)
  integer, intent(out) :: ierr

  real(dp), parameter :: eps_small = 1d-12  ! a very small number
  real(dp), parameter :: tol = 1d-8  ! fractional tolerance for bisection method
  real(dp) :: M2, M1dot, a, q, log_q
  real(dp) :: xL1, xL2, mu
  real(dp) :: Rd_over_a, PhiL1_dimless, PhiL2_dimless, PhiRd_dimless
  real(dp) :: GM2, Rd, Phi_units, PhiL1, PhiL2, PhiRd, omega_K
  real(dp) :: c1, c2, c3, c4

  integer :: i
  integer, parameter :: n_the = 50  ! number of grid points for disk thickness search
  real(dp), parameter :: the_grid_min = 0.1_dp
  real(dp), parameter :: the_grid_max = 1.0_dp
  real(dp) :: the_grid(n_the)  ! grid for disk thickness
  real(dp) :: T_arr(n_the)
  real(dp), parameter :: T_floor = 3.0d3  ! [K] the minimum value for disk temperature solution
  real(dp) :: T, T_max, the, the_left, the_right, the_min, the_max, separation_factor, dlogthe
  real(dp) :: T_left, T_right, f1, f1_left, f2, f2_left, f2_right, f, f_left
  real(dp) :: logT_arr(n_the), logthe_grid(n_the)

  ierr = 0

  ! key parameters
  M2 = accretor_mass * Msun
  M1dot = mass_transfer_rate * Msun/secyer
  a = orbital_separation * Rsun
  q = accretor_mass / donor_mass  ! mass ratio M2/M1

  log_q = log10(q)

  ! positions of Lagrangian points (based on analytic fits)
  xL1 = -0.0355_dp * log_q**2 + 0.251_dp * abs(log_q) + 0.500_dp  ! [a = SMA]
  xL2 = 0.0756_dp * log_q**2 - 0.424_dp * abs(log_q) + 1.699_dp  ! [a]

  if (log_q > 0.0_dp) then  ! m2 is more massive
     xL1 = 1.0_dp - xL1
     xL2 = 1.0_dp - xL2
  end if
  mu = q / (1.0_dp + q)

  ! outer disk radius
  Rd_over_a = pow4(1.0_dp - xL1) / mu
  ! relavent potential energies
  PhiL1_dimless = -((1.0_dp - mu)/abs(xL1) + mu/abs(1.0_dp - xL1) + 0.5_dp*(xL1 - mu)**2)  ! [G(M1+M2)/a]
  PhiL2_dimless = -((1.0_dp - mu)/abs(xL2) + mu/abs(1.0_dp - xL2) + 0.5_dp*(xL2 - mu)**2)  ! [G(M1+M2)/a]
  PhiRd_dimless = -(1.0_dp - mu + mu/Rd_over_a + 0.5_dp*(1.0_dp - mu)**2)

  GM2 = standard_cgrav * M2

  Rd = Rd_over_a*a
  Phi_units = standard_cgrav * (M2 / mu) / a
  PhiL1 = PhiL1_dimless * Phi_units
  PhiL2 = PhiL2_dimless * Phi_units
  PhiRd = PhiRd_dimless * Phi_units
  ! Keplerian frequency at Rd
  omega_K = sqrt(GM2 / pow3(Rd))


  ! constants involved in numerical solutions
  c1 = two_thirds * pi * crad * disk_alpha * Rd / (omega_K * M1dot)
  c2 = boltzm * Rd / (GM2 * disk_mu * mp)
  c3 = 8.0_dp * pi2 * crad * disk_alpha * clight * Rd**2 / (M1dot**2 * omega_K)
  c4 = 2.0_dp * pi * disk_mu * crad * disk_alpha * omega_K * mp * pow3(Rd) / (boltzm * M1dot)

  ! Create logarithmically spaced grid for grid search for disk thickness [only used at the beginning]
  do i = 1, n_the
     the_grid(i) = exp10(log10(the_grid_min) + (i - 1) * (log10(the_grid_max) - log10(the_grid_min)) / (n_the - 1))
  end do

  ! only T < T_max is possible to calculate
  T_max = pow(four_twentyseventh / (c1**2 * c2), one_ninth)

  do i = 1, n_the
     T_arr(i) = 0.0_dp
  end do

  do i = 1, n_the
     the = the_grid(i)
     ! use bisection method
     T_left = 0.1_dp * min(the**2 / c2, T_max)
     f1_left = f1_the_T_fL2(the, T_left, 0.0_dp, c1, c2)
     T_right = T_max
     do while (abs((T_left - T_right) / T_right) > tol)
        T = 0.5_dp * (T_left + T_right)
        f1 = f1_the_T_fL2(the, T, 0.0_dp, c1, c2)
        if (f1 * f1_left > 0) then
           T_left = T
           f1_left = f1
        else
           T_right = T
        end if
     end do
     T_arr(i) = 0.5_dp * (T_left + T_right)
  end do
  ! now we have obtained numerical relation between the and T
  do i = 1, n_the
     logT_arr(i) = log10(T_arr(i))
     logthe_grid(i) = log10(the_grid(i))
  end do
  dlogthe = logthe_grid(2) - logthe_grid(1)

  ! bisection to find the numerical solution to f2(the, T, fL2=0)=0
  the_right = 1.0_dp
  f2_right = f2_the_T_fL2(the_right, T_the_nofL2(the_right, logthe_grid, logT_arr, n_the, dlogthe), 0.0_dp, &
       c3, c4, M1dot, disk_alpha, omega_K, Rd, PhiRd, GM2, PhiL1, PhiL2)
  separation_factor = 0.95_dp
  the_left = separation_factor * the_right
  f2_left = f2_the_T_fL2(the_left, T_the_nofL2(the_left, logthe_grid, logT_arr, n_the, dlogthe), 0.0_dp, &
       c3, c4, M1dot, disk_alpha, omega_K, Rd, PhiRd, GM2, PhiL1, PhiL2)
  do while (f2_left * f2_right > 0)
     ! need to decrease the_left
     the_right = the_left
     f2_right = f2_left
     the_left = the_left * separation_factor
     f2_left = f2_the_T_fL2(the_left, T_the_nofL2(the_left, logthe_grid, logT_arr, n_the, dlogthe), 0.0_dp, &
          c3, c4, M1dot, disk_alpha, omega_K, Rd, PhiRd, GM2, PhiL1, PhiL2)
  end do
  ! now the solution is between the_left and the_right
  do while (abs((the_left - the_right) / the_right) > tol)
     the = 0.5_dp * (the_left + the_right)
     f2 = f2_the_T_fL2(the, T_the_nofL2(the, logthe_grid, logT_arr, n_the, dlogthe), 0.0_dp, &
          c3, c4, M1dot, disk_alpha, omega_K, Rd, PhiRd, GM2, PhiL1, PhiL2)
     if (f2 * f2_left > 0.0_dp) then
        the_left = the
        f2_left = f2
     else
        the_right = the
     end if
  end do
  ! solution
  the = 0.5_dp * (the_left + the_right)
  T = T_the_nofL2(the, logthe_grid, logT_arr, n_the, dlogthe)
  the_max = sqrt(three_eigths * c2 * T + one_fourth * (PhiL2 - PhiRd) / (GM2 / Rd) - one_eigth)

  if (the < the_max) then
     ! return a tiny numner
     fL2 = eps_small
  else
     the_min = ( 0.5_dp * sqrt((PhiL2 - PhiRd) / (GM2 / Rd) - 0.5_dp) )  ! corresponding to fL2=1, T=0
     ! need to find the maximum corresponding to fL2=0
     ! this is given by the intersection between T_the(the), T_the_nofL2(the)
     the_left = the_min
     the_right = 1.0_dp
     f_left = T_the(the_left, c2, PhiL2, PhiRd, GM2, Rd) - T_the_nofL2(the_left, logthe_grid, logT_arr, n_the, dlogthe)
     do while (abs((the_left - the_right) / the_right) > tol)
        the = 0.5_dp * (the_left + the_right)
        f = T_the(the, c2, PhiL2, PhiRd, GM2, Rd) - T_the_nofL2(the, logthe_grid, logT_arr, n_the, dlogthe)
        if (f * f_left > 0.0_dp) then
           the_left = the
           f_left = f
        else
           the_right = the
        end if
     end do
     the_max = 0.5_dp * (the_left + the_right)  ! this corresponds to fL2=0

     ! --- numerical solution for f2(the, T, fL2)=0 under non-zero fL2

     ! -- do not use exactly the_min (corresponding to T = 0, bc. kap table breaks down)
     ! -- define another the_min based on T_floor (kap table won't be a problem)
     the_min = sqrt( three_eigths * c2 * T_floor &
          + one_fourth * (PhiL2 - PhiRd) / (GM2 / Rd) - one_eigth )
     the_left = the_min
     f2_left = f2_the_T_fL2(the_left, &
          T_the(the_left, c2, PhiL2, PhiRd, GM2, Rd), &
          fL2_the(the_left, c1, c2, PhiL2, PhiRd, GM2, Rd), &
          c3, c4, M1dot, disk_alpha, omega_K, Rd, PhiRd, GM2, PhiL1, PhiL2)
     the_right = the_max / (1.0_dp + eps_small)
     ! bisection again
     do while (abs((the_left - the_right) / the_right) > tol)
        the = 0.5_dp * (the_left + the_right)
        f2 = f2_the_T_fL2(the, T_the(the, c2, PhiL2, PhiRd, GM2, Rd), fL2_the(the, c1, c2, PhiL2, PhiRd, GM2, Rd), &
             c3, c4, M1dot, disk_alpha, omega_K, Rd, PhiRd, GM2, PhiL1, PhiL2)
        if (f2 * f2_left > 0.0_dp) then
           the_left = the
           f2_left = f2
        else
           the_right = the
        end if
     end do
     ! solution
     the = 0.5_dp * (the_left + the_right)
     fL2 = fL2_the(the, c1, c2, PhiL2, PhiRd, GM2, Rd)
  end if

end subroutine eval_L2_mass_loss_fraction


! Helper Functions
real(dp) function f1_the_T_fL2(the, T, fL2, c1, c2)
  real(dp), intent(in) :: the, T, fL2, c1, c2
  f1_the_T_fL2 = c1 * pow4(T) * pow3(the) / (1.0_dp - fL2) - the**2 + c2 * T
end function f1_the_T_fL2

real(dp) function  kap(rho, T)
  ! simplified Kramers rule (cgs; approximate)
  real(dp), intent(in) :: rho, T
  kap = 0.34_dp + 3.0d24 * rho * pow(T, -3.5_dp)
end function kap

real(dp) function f2_the_T_fL2(the, T, fL2, &
     c3, c4, M1dot, disk_alpha, omega_K, Rd, PhiRd, GM2, PhiL1, PhiL2)
  real(dp), intent(in) :: the, T, fL2
  real(dp), intent(in) :: c3, c4, M1dot, disk_alpha, omega_K, Rd, PhiRd, GM2, PhiL1, PhiL2
  real(dp) :: x, U_over_P, rho

  real(dp), parameter :: seven_fourths = 7.0_dp / 4.0_dp

  x = c4 * pow3(T * the) / (1.0_dp - fL2)
  U_over_P = (1.5_dp + x) / (1.0_dp + one_third * x)
  rho = (1.0_dp - fL2) * M1dot / (2.0_dp * pi * disk_alpha * omega_K * pow3(Rd)) / pow2(the)
  f2_the_T_fL2 = &
       seven_fourths &
       - (1.5_dp * U_over_P + c3 * pow4(T) / kap(rho, T) / (1.0_dp - fL2) ** 2) * the**2 &
       - PhiRd / (GM2 / Rd) &
       + (PhiL1 - fL2 * PhiL2) / (GM2 / Rd) / (1.0_dp - fL2)
end function f2_the_T_fL2

real(dp) function T_the_nofL2(the, logthe_grid, logT_arr, n_the, dlogthe)
  real(dp), intent(in) :: the
  integer, intent(in) :: n_the
  real(dp), intent(in) :: logthe_grid(n_the), logT_arr(n_the)
  real(dp), intent(in) :: dlogthe
  real(dp) :: logthe, slope, logT
  integer :: i_the
  ! under the assumption fL2=0
  logthe = log10(the)
  if (logthe > logthe_grid(n_the-1)) then
     ! use analytic extrapolation
     T_the_nofL2 = exp10(logT_arr(n_the-1) - 0.25_dp * (logthe - logthe_grid(n_the-1)))
  else if (logthe < logthe_grid(1)) then
     ! analytic extrapolation
     T_the_nofL2 = exp10(logT_arr(1) + 2.0_dp * (logthe - logthe_grid(1)))
  else
     i_the = floor((logthe - logthe_grid(1)) / dlogthe) + 1
     slope = (logT_arr(i_the + 1) - logT_arr(i_the)) / dlogthe
     logT = logT_arr(i_the) + (logthe - logthe_grid(i_the)) * slope
     T_the_nofL2 = exp10(logT)
  end if
end function T_the_nofL2

real(dp) function T_the(the, c2, PhiL2, PhiRd, GM2, Rd)
  real(dp), intent(in) :: the
  real(dp), intent(in) :: c2, PhiL2, PhiRd, GM2, Rd
  ! only for non-zero fL2
  T_the = (8.0_dp * the**2 + 1.0_dp - 2.0_dp * (PhiL2 - PhiRd) / (GM2 / Rd)) / (3.0_dp * c2)
end function T_the


real(dp) function fL2_the(the, c1, c2, PhiL2, PhiRd, GM2, Rd)
  real(dp), intent(in) :: the
  real(dp), intent(in) :: c1, c2, PhiL2, PhiRd, GM2, Rd
  real(dp) :: T
  ! only for non-zero fL2
  T = T_the(the, c2, PhiL2, PhiRd, GM2, Rd)
  fL2_the =  1.0_dp - c1 * pow4(T) * pow3(the) / (the**2 - c2 * T)
end function fL2_the

! ----------------------------------------------------------------------------------------------
! helper functions

real(dp) function get_min_r(binary_id) result(min_r)
  integer, intent(in) :: binary_id
  integer :: ierr
  type(binary_info), pointer :: b
  real(dp) :: qratio
  ierr = 0
  min_r = 0
  call binary_ptr(binary_id, b, ierr)
  if (ierr /= 0) then
     write(*,*) 'failed in binary_ptr'
     STOP "get_min_r can't get binary_pointer"
  end if

  ! reimplement Ulrich & Burger 1976 's fit to Lubow & Shu 1975
  ! as in de Mink+13, copying code from $MESA_DIR/binary/private/binary_mdot.f90
  qratio = b% m(b% a_i) / b% m(b% d_i)
  qratio = min(max(qratio,0.0667d0),15d0)
  min_r = 0.0425d0*b% separation*pow(qratio+qratio*qratio, 0.25d0) ! same units as separation, cm
end function get_min_r

! this is copied from $MESA_DIR/binary/private/binary_wind.f90
subroutine eval_wind_xfer_fractions(binary_id, ierr)
  integer, intent(in) :: binary_id
  integer, intent(out) :: ierr
  type (binary_info), pointer :: b
  ierr = 0
  call binary_ptr(binary_id, b, ierr)
  if (ierr /= 0) then
     write(*,*) 'failed in binary_ptr'
     return
  end if

  ! for the primary
  if (b% point_mass_i /= 1) then
     if (.not. b% do_wind_mass_transfer_1 .or. b% model_twins_flag) then
        b% wind_xfer_fraction(1) = 0d0
     else if(.not. b% use_other_binary_wind_transfer) then
        call Bondi_Hoyle_wind_transfer(b% binary_id, 1, ierr)
        if (ierr /=0) then
           write(*,*) "Error in Bondi_Hoyle_wind_transfer(b% binary_id, 1, ierr)"
           return
        end if
     else
        call b% other_binary_wind_transfer(b% binary_id, 1, ierr)
        if (ierr /=0) then
           write(*,*) "Error in other_binary_wind_transfer(b% binary_id, 1, ierr)"
           return
        end if
     end if
  end if

  ! check if secondary needs wind transfer
  if (b% point_mass_i /= 2) then
     if (.not. b% do_wind_mass_transfer_2) then
        b% wind_xfer_fraction(2) = 0d0
     else if(.not. b% use_other_binary_wind_transfer) then
        call Bondi_Hoyle_wind_transfer(b% binary_id, 2, ierr)
        if (ierr /=0) then
           write(*,*) "Error in Bondi_Hoyle_wind_transfer(b% binary_id, 2, ierr)"
           return
        end if
     else
        call b% other_binary_wind_transfer(b% binary_id, 2, ierr)
        if (ierr /=0) then
           write(*,*) "Error in other_binary_wind_transfer(b% binary_id, 2, ierr)"
           return
        end if
     end if
  end if

end subroutine eval_wind_xfer_fractions

! this is copied from $MESA_DIR/binary/private/binary_wind.f90
subroutine Bondi_Hoyle_wind_transfer(binary_id, s_i, ierr)
  integer, intent(in) :: binary_id, s_i ! s_i is index of the wind mass losing star
  integer, intent(out) :: ierr

  ! wind transfer fraction based on Bondi-Hoyle mechanism as described in
  ! Hurley et al. 2002, MNRAS, 329, 897-928

  type(binary_info), pointer :: b
  type (star_info), pointer :: s
  real(dp) :: v_orb, v_wind
  real(dp) :: alpha, beta  ! Bondi-Hoyle alpha, beta for that star
  real(dp) :: max_xfer  ! Maximum transfer fraction

  call binary_ptr(binary_id, b, ierr)
  if (ierr /= 0) then
     write(*,*) 'failed in binary_ptr'
     return
  end if

  if (s_i == 1) then
     s => b% s1
     alpha = b% wind_BH_alpha_1
     beta = b% wind_BH_beta_1
     max_xfer = b% max_wind_transfer_fraction_1
  else
     s => b% s2
     alpha = b% wind_BH_alpha_2
     beta = b% wind_BH_beta_2
     max_xfer = b% max_wind_transfer_fraction_2
  end if

  ! orbital speed Hurley et al 2002 eq. 8
  v_orb = sqrt(standard_cgrav * (b% m(1) + b% m(2)) / b% separation) !cm/s

  ! windspeed from Hurley et al 2002 eq. 9
  v_wind = sqrt(2d0 * beta *  standard_cgrav * b% m(s_i) / b% r(s_i))

  ! Bondi-Hoyle transfer fraction Hurley et al. 2002 eq. 6
  b% wind_xfer_fraction(s_i) = alpha / pow2(b% separation) /&
       (2d0 * sqrt(1d0 - pow2(b% eccentricity))) *&
       pow2(standard_cgrav * b% m(3-s_i) / pow2(v_wind)) *&
       pow(1d0 + pow2(v_orb/v_wind),-1.5d0)

  ! limit to provided maximum
  b% wind_xfer_fraction(s_i) = min(max_xfer, b% wind_xfer_fraction(s_i))

end subroutine Bondi_Hoyle_wind_transfer
! ----------------------------------------------------------------------------------------------
! The magic happens down here

! other_accreted_material_j
! This routine implements a simplified version of the solution
! for the boundary layer between an accretor star and an accretion disk
! following Paczsynski 1991. As the accretor approaches critical rotation
! the amount of accreted specific angular momentum (w.r.t. the center of
! the accretor) goes from the Keplerian surface value sqrt(GM2R2) to zero.
! For rotation approaching critical a smooth tanh interpolation is
! used between the two values
subroutine disk_accreted_material_j(binary_id, ierr)
  use binary_def, only : binary_info, binary_ptr
  use const_def, only: dp
  use utils_lib, only: mesa_error
  integer, intent(in) :: binary_id
  integer, intent(out) :: ierr
  real(dp) :: min_r, omega_ratio_acc
  real(dp) :: j_acc_disk, j_acc_crit, j_interp
  type (binary_info), pointer :: b
  ierr = 0
  call binary_ptr(binary_id, b, ierr)
  if (ierr /= 0) then
     write(*,*) 'failed in binary_ptr'
     return
  end if

  ! use disk AM until a certain omega/omega_crit, then ramp up to zero quickly
  ! this mimics qualitatively Paczinski 1991's solution
  j_acc_disk = sqrt(standard_cgrav * b% m(b% a_i) * b% r(b% a_i)) ! specific Keplerian AM at accretor surface w.r.t. accretor
  j_acc_crit = 0.0d0 ! no AM accretion close to critical rotation, as per Paczynski's work

  omega_ratio_acc = b% s_accretor% omega_avg_surf / b% s_accretor% omega_crit_avg_surf
  print *, "omega_div_omega_crit accretor = ", omega_ratio_acc, "j_keplerian2 = ", j_acc_disk

  if (b% accretion_mode == 1) then
     ! direct impact
     min_r = get_min_r(binary_id)
     b% s_accretor% accreted_material_j = &
          sqrt(standard_cgrav * b% m(b% a_i) * 1.7d0*min_r)
  else if (b% accretion_mode == 2) then
     ! disk is present
     ! interpolate between Keplerian and 0.0d0 using omega/omega_crit as interpolation variable
     if (omega_ratio_acc <= b% x_ctrl(1)) then
        j_interp = j_acc_disk  ! Keplerian value
     else if (omega_ratio_acc >= b% x_ctrl(2)) then
        j_interp = j_acc_crit  ! null -- or implement a solution of Paczynski's idea with dependence on b% mtransfer_rate of b% s_accretor% mstar_dot
     else
        j_interp = j_acc_disk + &
             0.5d0*(j_acc_crit-j_acc_disk) * &
             (1+tanh(b%x_ctrl(3)* &  ! tanh for smooth interpolation
             (omega_ratio_acc - ((b%x_ctrl(1) + b%x_ctrl(2))/2.0) / &
             (b%x_ctrl(2) - b%x_ctrl(1)))))
     end if
     b% s_accretor% accreted_material_j = j_interp
  else if (b% accretion_mode == 3) then
     ! implement what jdot to use for direct impact at critical rotation
     b% s_accretor% accreted_material_j = j_acc_crit ! assume no AM accretion
  end if
  ! save for reuse later
  ! b% xtra(1) = j_acc_disk
  ! b% xtra(2) = b% s_accretor% accreted_material_j
  ! print *, "j_accreted_2 - j_keplerian2", b% s_accretor% accreted_material_j - j_acc_disk
  ! fix ratio to Keplerian
  b% acc_am_div_kep_am = b% s_accretor% accreted_material_j / &
          sqrt(standard_cgrav * b% m(b% a_i) * b% r(b% a_i))
end subroutine disk_accreted_material_j


! other_adjust_mdots
! This routine implements the mass lost from the accretion disk
! following Lu et al. 2023 Based on an energetic argument, we
! calculate the fraction of the mass lost by the donor that escapes
! through L2 and implement both the decrease in accretion rate of the
! secondary and (with other_jdot) the loss of orbital AM the change in
! mass transfer due to the L2 losses are added as a term to the mass
! tranfer efficiency, which becomes 1-alpha-beta-delta-fL2 (see
! comment in $MESA_DIR/binary/defaults/controls.default)
subroutine L2_adjust_mdots(binary_id, ierr)
  use binary_def, only : binary_info, binary_ptr
  use binary_lib, only : binary_eval_accreted_material_j
  ! use binary_wind, only: eval_wind_xfer_fractions
  use const_def, only: dp,  standard_cgrav
  real(dp) :: qratio, min_r ! needed for reimplementation of UB76 fit to LS75
  integer, intent(in) :: binary_id
  integer, intent(out) :: ierr
  real(dp) :: fL2, disk_alpha, actual_mtransfer_rate
  type (binary_info), pointer :: b
  ierr = 0
  call binary_ptr(binary_id, b, ierr)
  if (ierr /= 0) then
     write(*,*) 'failed in binary_ptr'
     return
  end if

  ! initialize
  actual_mtransfer_rate = 0d0
  ! fL2 = 0.0d0

  b% fixed_xfer_fraction = 1 - b% mass_transfer_alpha - b% mass_transfer_beta - &
       b% mass_transfer_delta

  ! turn off Eddington limit for non-compact objects
  b% mdot_edd = 0.0d0
  b% mdot_edd_eta = 0d0
  ! solve wind mass transfer
  ! b% mdot_wind_transfer(b% d_i) is a negative number that gives the
  ! amount of mass transferred by unit time from the donor to the
  ! accretor.
  call eval_wind_xfer_fractions(b% binary_id, ierr)
  if (ierr/=0) then
     write(*,*) "Error in eval_wind_xfer_fractions"
     return
  end if

  ! remove wind fractions
  b% mdot_wind_transfer(b% d_i) = b% s_donor% mstar_dot * &
       b% wind_xfer_fraction(b% d_i)
  b% mdot_wind_transfer(b% a_i) = b% s_accretor% mstar_dot * &
       b% wind_xfer_fraction(b% a_i)


  ! -----------------------------------------------------------------------
  ! Set mdot for the donor
  b% s_donor% mstar_dot = b% s_donor% mstar_dot + b% mtransfer_rate - &
       b% mdot_wind_transfer(b% a_i) ! accretion from the wind of the accretor
  ! -----------------------------------------------------------------------

  ! -----------------------------------------------------------------------
  ! Set mdot for accretor
  ! First, check for disk
  min_r = get_min_r(binary_id)
  ! print *, "min_r", min_r, "r_acc", b% r(b% a_i), "disk?", b% r(b% a_i) < min_r
  if (b% r(b% a_i) < min_r) then
     b% accretion_mode = 2 ! means there is a disk
     ! print *, "There is a disk!"
     disk_alpha = 0.1d0
     call eval_L2_mass_loss_fraction(b% s_donor% m(1)/Msun, &
          b% s_accretor% m(1)/Msun, &
          b% mtransfer_rate/(Msun/secyer), &
          b% separation/Rsun, &
          disk_alpha, &
          b% s_donor% mu(1), &
          fL2, ierr)
  else
     ! print *, "Direct impact!"
     if (b% s_accretor% omega_avg_surf / b% s_accretor% omega_crit_avg_surf >= b% x_ctrl(1)) then
        b% accretion_mode = 3 ! 3 means there is a decretion disk
        ! print *, "Direct impact and a decretion disk?!", b% s_accretor% omega_avg_surf / b% s_accretor% omega_crit_avg_surf
     else
        ! direct impact, no disk, no update to fL2
        b% accretion_mode = 1
     end if
  end if

  ! turn off accretion luminosity for non-compact stars
  b% accretion_luminosity = 0.0d0

  ! modify L2 mass loss value:
  b% s_accretor% mstar_dot = b% s_accretor% mstar_dot - &
       b% mtransfer_rate*(b% fixed_xfer_fraction - fL2) - b% mdot_wind_transfer(b% d_i)
  ! -----------------------------------------------------------------------

  ! now deal with angular momentum
  if (b% use_other_accreted_material_j .eqv. .true.) then
     call b% other_accreted_material_j(b% binary_id, ierr)
  else
     call binary_eval_accreted_material_j(b% binary_id, ierr)
  end if

  ! Finally adjust all mdots

  ! mdot_system_transfer is mass lost from the vicinity of each star
  ! due to inefficient rlof mass transfer, mdot_system_cct is mass lost
  ! from a circumbinary coplanar toroid.
  b% mdot_system_transfer(b% d_i) = b% mtransfer_rate * b% mass_transfer_alpha

  ! combine RLOF and winds
  actual_mtransfer_rate = b% mtransfer_rate * (b% fixed_xfer_fraction - fL2) + &
       b% mdot_wind_transfer(b% d_i) !defined negative
  b% component_mdot(b% a_i) = -actual_mtransfer_rate

  b% mdot_system_cct = b% mtransfer_rate * b% mass_transfer_delta

  b% mdot_system_transfer(b% a_i) = (actual_mtransfer_rate + b% component_mdot(b% a_i)) &
       + b% mtransfer_rate * b% mass_transfer_beta
end subroutine L2_adjust_mdots

! other_extra_jdot
subroutine L2_extra_jdot(binary_id, ierr)
  use binary_def, only : binary_info, binary_ptr
  use const_def, only: dp, pi
  integer, intent(in) :: binary_id
  integer, intent(out) :: ierr
  real(dp) :: fL2, disk_alpha, xL2, log_q, omega_orb
  real(dp) :: L2_am_loss_rate !, delta_j, disk_torque_on_orbit
  real(dp), parameter  :: g_b = 1.0d0
  type (binary_info), pointer :: b
  ierr = 0
  call binary_ptr(binary_id, b, ierr)
  if (ierr /= 0) then
     write(*,*) 'failed in binary_ptr'
     return
  end if
  ! initialize
  b% extra_jdot = 0.0d0
  ! fL2 = 0.0d0
  L2_am_loss_rate = 0.0d0
  ! disk_torque_on_orbit = 0.0d0

  if ((b% accretion_mode == 1)) then
     ! no disk, no extra torque
     return
  else if (b% accretion_mode == 2) then
     ! disk (accretion if 2, decretion+direct impact if 3)
     disk_alpha = 0.1d0
     call eval_L2_mass_loss_fraction(b% s_donor% m(1)/Msun, &
          b% s_accretor% m(1)/Msun, &
          b% mtransfer_rate/(Msun/secyer), &
          b% separation/Rsun, &
          disk_alpha, &
          b% s_donor% mu(1), &
          fL2, ierr)
     if (fL2 == 0.0d0) then
        ! there is a disk but no L2 mass loss
        return
     else
        ! there is a disk and there is L2 mass loss
        log_q = log10(b% s_accretor% m(1)/ b% s_donor% m(1))
        ! position of L2 w.r.t. center of mass according to Lu et al. 23 fit in units of separation
        xL2 = 0.0756_dp * log_q**2 - 0.424_dp * abs(log_q) + 1.699_dp  ! xL2 = rL2/a
        omega_orb = 2*pi/b% period ! 1/sec
        ! add AM loss from L2, see Eq. 36 in Lu et al. 2023
        L2_am_loss_rate = ((b% s_donor%mstar_dot * fL2)) * g_b * & ! amount of mass lost at L2, note that mstar_dot is negative
         omega_orb * pow2(b% separation) * pow2(xL2-(b% s_accretor% m(1) / (b%s_accretor% m(1)+ b% s_donor% m(1))))
        ! print *, "L2 am loss", L2_am_loss_rate
        ! delta_j = b% xtra(2) - b% xtra(1) ! accreted by secondary minus Keplerian at secondary surface
        !                                   ! if negative, assume it goes in the secondary accretion disk, and feed it
        !                                   ! back to the orbit
        ! print *, "delta_j", delta_j
        ! if (delta_j < 0.0) then
        !    ! accretor accreted less less AM than Keplerian
        !    b% extra_jdot = L2_am_loss_rate - delta_j/b% s_donor% dt ! delta_j is negative, -delta_j/timestep is a positive torque
        ! else
        !    ! accretor got >= keplerian, nothing to add
        !    b% extra_jdot = L2_am_loss_rate
        ! end if
        b% extra_jdot = L2_am_loss_rate
     end if
  else if (b% accretion_mode == 3) then
       ! decretion disk and direct hit, treat decretion disk like accretion disk
     return
  end if

  print *, "extra_jdot = ", b% extra_jdot
end subroutine L2_extra_jdot
