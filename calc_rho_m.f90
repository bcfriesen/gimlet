! Calculate mass density from baryon density and dark matter density.

! NOTE: this routine reproduces the results of old gimlet byte for byte,
! but I think it's wrong. In particular, it adds together the baryon and
! dark matter densities, with each in their own mean units, i.e., it
! adds together terms like (rho_b / rho_b_bar) + (rho_dm / rho_dm_bar),
! where rho_b_bar /= rho_dm_bar.

subroutine calc_rho_m(rho_b, rho_dm, lo, hi, rho_b_ng, rho_dm_ng, rho_m_ng, &
                      omega_b, omega_m, mean_density, mean_dm_density, rho_m)
  use, intrinsic :: iso_c_binding
  implicit none

  integer(c_int), intent(in) :: lo(3), hi(3), rho_b_ng, rho_dm_ng, rho_m_ng
  real(c_double), intent(in) :: rho_b  (lo(1)-rho_b_ng:hi(1)+rho_b_ng, &
                                        lo(2)-rho_b_ng:hi(2)+rho_b_ng, &
                                        lo(3)-rho_b_ng:hi(3)+rho_b_ng, 1)
  real(c_double), intent(in) :: rho_dm (lo(1)-rho_dm_ng:hi(1)+rho_dm_ng, &
                                        lo(2)-rho_dm_ng:hi(2)+rho_dm_ng, &
                                        lo(3)-rho_dm_ng:hi(3)+rho_dm_ng, 1)
  real(c_double), intent(in) :: omega_b, omega_m, mean_density, mean_dm_density

  real(c_double), intent(out) :: rho_m (lo(1)-rho_m_ng:hi(1)+rho_m_ng, &
                                        lo(2)-rho_m_ng:hi(2)+rho_m_ng, &
                                        lo(3)-rho_m_ng:hi(3)+rho_m_ng, 1)

  real(c_double) :: baryon_fac, dm_fac
  integer(c_int) :: i, j, k

  baryon_fac = omega_b / omega_m
  dm_fac = 1.0 - baryon_fac

  do k = lo(3), hi(3)
    do j = lo(2), hi(2)
      do i = lo(1), hi(1)
        rho_m (i, j, k, 1) = baryon_fac * rho_b (i, j, k, 1) / mean_density + &
                             dm_fac * rho_dm(i, j, k, 1) / mean_dm_density
      end do
    end do
  end do

end subroutine calc_rho_m
