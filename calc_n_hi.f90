! Calculate H I density using the Nyx EOS.

subroutine calc_n_hi(z, density, e_int, mean_density, omega_b, h_0, lo, hi, state_ng, n_hi_ng, n_hi)
  use, intrinsic :: iso_c_binding
  use eos_module, only: nyx_eos_nh0_and_nhep
  use fundamental_constants_module, only: e_to_cgs
  implicit none

  integer(c_int), intent(in) :: lo(3), hi(3), state_ng, n_hi_ng
  real(c_double), intent(in) :: z, mean_density, omega_b, h_0
  real(c_double), intent(in) :: density (lo(1)-state_ng:hi(1)+state_ng, &
                                         lo(2)-state_ng:hi(2)+state_ng, &
                                         lo(3)-state_ng:hi(3)+state_ng)
  real(c_double), intent(in) :: e_int (lo(1)-state_ng:hi(1)+state_ng, &
                                       lo(2)-state_ng:hi(2)+state_ng, &
                                       lo(3)-state_ng:hi(3)+state_ng)

  real(c_double), intent(out) :: n_hi (lo(1)-n_hi_ng:hi(1)+n_hi_ng, &
                                       lo(2)-n_hi_ng:hi(2)+n_hi_ng, &
                                       lo(3)-n_hi_ng:hi(3)+n_hi_ng)

  real(c_double), parameter :: rho_crit_100_cgs = 1.8788200387d-29

  real(c_double) :: comoving_a, a3, h
  real(c_double) :: cell_n_hi, cell_nhep, e_int_cell, rho_cell, rho_eos_units
  integer :: i, j, k

  h = h_0 / 100.0 ! Hubble parameter
  comoving_a = 1.0 / (1.0 + z)
  a3 = comoving_a*comoving_a*comoving_a

  do k = lo(3), hi(3)
    do j = lo(2), hi(2)
      do i = lo(1), hi(1)

        rho_cell = density (i, j, k)
        e_int_cell = e_int (i, j, k) / rho_cell

        rho_eos_units = rho_cell * omega_b * h*h * rho_crit_100_cgs / (a3 * mean_density)

        call nyx_eos_nh0_and_nhep(z, rho_eos_units, e_int_cell*e_to_cgs, cell_n_hi, cell_nhep)

        n_hi (i, j, k) = cell_n_hi

      end do
    end do
  end do
end subroutine calc_n_hi
