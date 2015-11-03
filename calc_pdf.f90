! Calculate PDF of a scalar field.

subroutine calc_pdf(dist_mf, lo, hi, ng, num_bins, bin_edges, bin_count, bin_x_sum)
    use, intrinsic :: iso_c_binding
    implicit none

    integer(c_int), intent(in) :: lo(3), hi(3), ng
    real(c_double), intent(in) :: dist_mf (lo(1)-ng:hi(1)+ng, lo(2)-ng:hi(2)+ng, lo(3)-ng:hi(3)+ng, 1)
    integer(c_int), intent(in) :: num_bins
    real(c_double), intent(in) :: bin_edges (num_bins+1)

    integer(c_int), intent(out) :: bin_count (num_bins)
    real(c_double), intent(out) :: bin_x_sum (num_bins)

    real(c_double) :: x_min, x_max
    real(c_double) :: xi

    integer(c_int) :: i, j, k
    integer(c_int) :: bin_index, jj

    x_min = minval(bin_edges)
    x_max = maxval(bin_edges)

    bin_count (:) = 0
    bin_x_sum (:) = 0.0d0

    do k = lo(3), hi(3)
      do j = lo(2), hi(2)
        do i = lo(1), hi(1)

          xi = dist_mf(i, j, k, 1)
          if (xi < x_min .or. xi > x_max) cycle

          ! Find the bin each value belongs in.
          bin_index = -1
          do jj = 1, num_bins
            if (xi >= bin_edges(jj) .and. xi < bin_edges(jj+1)) then
                bin_index = jj
                exit
            end if
          end do
          ! If value exactly equals the largest bin edge, put it in the largest bin.
          if (bin_index == -1) bin_index = num_bins

          bin_count (bin_index) = bin_count (bin_index) + 1
          bin_x_sum (bin_index) = bin_x_sum (bin_index) + xi

        end do
      end do
    end do

end subroutine calc_pdf
