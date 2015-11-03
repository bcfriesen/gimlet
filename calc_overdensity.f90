! Given a scalar field f(x), calculate g(x) = (x/x_bar) - 1, where x_bar is the mean. "Overdensity" is a misnomer because this can
! be used for any scalar field, not just density. But I don't know what the technical term for this is.

subroutine calc_overdensity(density, lo, hi, ng, mean_density_inv, overdensity)
    use, intrinsic :: iso_c_binding
    implicit none

    integer(c_int), intent(in) :: lo(3), hi(3), ng
    real(c_double), intent(in) :: density (lo(1)-ng:hi(1)+ng, lo(2)-ng:hi(2)+ng, lo(3)-ng:hi(3)+ng, 1)
    real(c_double), intent(in) :: mean_density_inv

    real(c_double), intent(out) :: overdensity (lo(1)-ng:hi(1)+ng, lo(2)-ng:hi(2)+ng, lo(3)-ng:hi(3)+ng, 1)

    integer :: i, j, k

    do k = lo(3), hi(3)
      do j = lo(2), hi(2)
        do i = lo(1), hi(1)
          overdensity (i, j, k, 1) = density (i, j, k, 1) * mean_density_inv - 1.0
        end do
      end do
    end do

end subroutine calc_overdensity
