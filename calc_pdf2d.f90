! Calculate PDF of two scalar fields along two axes.

subroutine calc_pdf2d(field1, field2, lo, hi, ng, &
                      num_x_bins, num_y_bins, bin_edges_x, bin_edges_y, bin_count, &
                      bin_x_sum, bin_y_sum)
    use, intrinsic :: iso_c_binding
    implicit none

    ! Assumes the two MultiFabs have the same shape and # of ghost cells.

    integer(c_int), intent(in) :: lo(3), hi(3), ng
    real(c_double), intent(in) :: field1 (lo(1)-ng:hi(1)+ng, lo(2)-ng:hi(2)+ng, lo(3)-ng:hi(3)+ng, 1)
    real(c_double), intent(in) :: field2 (lo(1)-ng:hi(1)+ng, lo(2)-ng:hi(2)+ng, lo(3)-ng:hi(3)+ng, 1)
    integer(c_int), intent(in) :: num_x_bins, num_y_bins
    real(c_double), intent(in) :: bin_edges_x (num_x_bins+1), bin_edges_y (num_y_bins+1)

    integer(c_int), intent(out) :: bin_count (num_x_bins*num_y_bins)
    real(c_double), intent(out) :: bin_x_sum (num_x_bins*num_y_bins), bin_y_sum (num_x_bins*num_y_bins)

    real(c_double) :: x_min, x_max, y_min, y_max
    real(c_double) :: xi, yi

    integer(c_int) :: i, j, k
    integer(c_int) :: bin_index_x, bin_index_y, jj
    integer(c_int) :: num_bins
    integer(c_int) :: ibin

    num_bins = num_x_bins*num_y_bins

    x_min = minval(bin_edges_x)
    x_max = maxval(bin_edges_x)
    y_min = minval(bin_edges_y)
    y_max = maxval(bin_edges_y)

    bin_count (:) = 0
    bin_x_sum (:) = 0.0d0
    bin_y_sum (:) = 0.0d0

    do k = lo(3), hi(3)
      do j = lo(2), hi(2)
        do i = lo(1), hi(1)

          xi = field1(i, j, k, 1)
          yi = field2(i, j, k, 1)

          if (xi < x_min .or. xi > x_max .or. yi < y_min .or. yi > y_max) cycle

          ! Dind the x bin.
          bin_index_x = -1
          do jj = 1, num_x_bins
            if (xi >= bin_edges_x(jj) .and. xi < bin_edges_x(jj+1)) then
                bin_index_x = jj
                exit
            end if
          end do
          ! If we can't find the bin, assume it goes in the last one.
          if (bin_index_x == -1) bin_index_x = num_x_bins

          ! Find the y bin.
          bin_index_y = -1
          do jj = 1, num_y_bins
            if (yi >= bin_edges_y(jj) .and. yi < bin_edges_y(jj+1)) then
                bin_index_y = jj
                exit
            end if
          end do
          ! If we can't find the bin, assume it goes in the last one.
          if (bin_index_y == -1) bin_index_y = num_y_bins

          ! Flatten index into (x, y) array.
          ibin = bin_index_y + ((bin_index_x-1) * num_y_bins)

          bin_count (ibin) = bin_count (ibin) + 1
          bin_x_sum (ibin) = bin_x_sum (ibin) + xi
          bin_y_sum (ibin) = bin_y_sum (ibin) + yi

        end do
      end do
    end do

end subroutine calc_pdf2d
