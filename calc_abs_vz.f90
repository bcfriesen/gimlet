! Calculate magnitude of z-component of velocity. NOTE the units are km/s, not cm/s.

subroutine calc_abs_vz(zmom, density, zmom_ng, density_ng, vz_ng, lo, hi, abs_vz)
  use, intrinsic :: iso_c_binding
  implicit none

  integer(c_int), intent(in) :: lo(3), hi(3), zmom_ng, density_ng, vz_ng
  real(c_double), intent(in) :: zmom (lo(1)-zmom_ng:hi(1)+zmom_ng, &
                                      lo(2)-zmom_ng:hi(2)+zmom_ng, &
                                      lo(3)-zmom_ng:hi(3)+zmom_ng, 1)
  real(c_double), intent(in) :: density (lo(1)-density_ng:hi(1)+density_ng, &
                                         lo(2)-density_ng:hi(2)+density_ng, &
                                         lo(3)-density_ng:hi(3)+density_ng, 1)

  real(c_double), intent(out) :: abs_vz (lo(1)-vz_ng:hi(1)+vz_ng, &
                                         lo(2)-vz_ng:hi(2)+vz_ng, &
                                         lo(3)-vz_ng:hi(3)+vz_ng, 1)

  integer :: i, j, k

  do k = lo(3), hi(3)
    do j = lo(2), hi(2)
      do i = lo(1), hi(1)
        abs_vz(i, j, k, 1) = abs(zmom(i, j, k, 1)) / density(i, j, k, 1)
      end do
    end do
  end do
end subroutine calc_abs_vz
