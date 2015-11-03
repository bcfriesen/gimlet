! Calculate flux from optical depth using the formula F ~ exp(-tau).

subroutine calc_flux(tau, lo, hi, ng, A, flux)
    use, intrinsic :: iso_c_binding
    implicit none

    integer(c_int), intent(in) :: lo(3), hi(3), ng
    real(c_double), intent(in) :: tau (lo(1)-ng:hi(1)+ng, lo(2)-ng:hi(2)+ng, lo(3)-ng:hi(3)+ng, 1)
    real(c_double), intent(in) :: A

    real(c_double), intent(out) :: flux (lo(1)-ng:hi(1)+ng, lo(2)-ng:hi(2)+ng, lo(3)-ng:hi(3)+ng, 1)

    integer(c_int) :: i, j, k

    do k = lo(3), hi(3)
      do j = lo(2), hi(2)
        do i = lo(1), hi(1)
          flux (i, j, k, 1) = exp (-A * tau (i, j, k, 1))
        end do
      end do
    end do

end subroutine calc_flux
