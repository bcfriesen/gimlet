! Calculate velocity magnitude from momenta and density. NOTE the resulting units are km/s, NOT cm/s.

subroutine calc_abs_v(xmom, ymom, zmom, mom_ng, density, density_ng, lo, hi, abs_v)
  use, intrinsic :: iso_c_binding
  implicit none

  integer(c_int), intent(in) :: lo(3), hi(3), mom_ng, density_ng
  real(c_double), intent(in) :: xmom (lo(1)-mom_ng:hi(1)+mom_ng, &
                                      lo(2)-mom_ng:hi(2)+mom_ng, &
                                      lo(3)-mom_ng:hi(3)+mom_ng, 1)
  real(c_double), intent(in) :: ymom (lo(1)-mom_ng:hi(1)+mom_ng, &
                                      lo(2)-mom_ng:hi(2)+mom_ng, &
                                      lo(3)-mom_ng:hi(3)+mom_ng, 1)
  real(c_double), intent(in) :: zmom (lo(1)-mom_ng:hi(1)+mom_ng, &
                                      lo(2)-mom_ng:hi(2)+mom_ng, &
                                      lo(3)-mom_ng:hi(3)+mom_ng, 1)
  real(c_double), intent(in) :: density (lo(1)-density_ng:hi(1)+density_ng, &
                                         lo(2)-density_ng:hi(2)+density_ng, &
                                         lo(3)-density_ng:hi(3)+density_ng, 1)

  real(c_double), intent(out) :: abs_v(lo(1)-mom_ng:hi(1)+mom_ng, &
                                       lo(2)-mom_ng:hi(2)+mom_ng, &
                                       lo(3)-mom_ng:hi(3)+mom_ng, 1)

  integer :: i, j, k
  real(c_double) :: xvel, yvel, zvel

  do k = lo(3), hi(3)
    do j = lo(2), hi(2)
      do i = lo(1), hi(1)
        xvel = xmom(i, j, k, 1) / density(i, j, k, 1)
        yvel = ymom(i, j, k, 1) / density(i, j, k, 1)
        zvel = zmom(i, j, k, 1) / density(i, j, k, 1)

        ! NOTE: THE UNITS ARE KM/S, NOT CM/S
        abs_v(i, j, k, 1) = sqrt(xvel*xvel + yvel*yvel + zvel*zvel)
      end do
    end do
  end do
end subroutine calc_abs_v
