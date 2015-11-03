! Calculate power spectrum P(k) of the 3-D Fourier transform of a scalar field (as opposed to the FT along 1-D pencils).

subroutine calc_ps3d (overdensity_fft_real, overdensity_fft_imag, lo, hi, ng, num_bins, k_bin_edges, domain_length, &
                      domain_grid_length, k_bin_count, k_bin_power_weighted_k_sum, k_bin_power_sum)

  use, intrinsic :: iso_c_binding
  implicit none

  integer(c_int), intent(in) :: lo(3), hi(3), ng, num_bins
  real(c_double), intent(in), target :: overdensity_fft_real (lo(1)-ng:hi(1)+ng, lo(2)-ng:hi(2)+ng, lo(3)-ng:hi(3)+ng, 1)
  real(c_double), intent(in), target :: overdensity_fft_imag (lo(1)-ng:hi(1)+ng, lo(2)-ng:hi(2)+ng, lo(3)-ng:hi(3)+ng, 1)
  real(c_double), intent(in) :: k_bin_edges (num_bins+1)
  real(c_double), intent(in) :: domain_length
  integer(c_int), intent(in) :: domain_grid_length

  integer(c_int), intent(out) :: k_bin_count (num_bins)
  real(c_double), intent(out) :: k_bin_power_weighted_k_sum (num_bins)
  real(c_double), intent(out) :: k_bin_power_sum (num_bins)

  integer(c_int) :: i, j, k
  integer(c_int) :: kk
  integer(c_int) :: k_bin_index
  integer(c_int) :: k_index_i, k_index_j, k_index_k

  real(c_double) :: k_fund
  real(c_double) :: k_mag
  real(c_double) :: k_max
  real(c_double) :: k_min
  real(c_double) :: p_i
  real(c_double) :: pi

  pi = acos(-1.0)

  k_fund = 2.0 * pi / domain_length

  k_min = minval(k_bin_edges)
  k_max = maxval(k_bin_edges)

  k_bin_count (:) = 0
  k_bin_power_weighted_k_sum (:) = 0.0
  k_bin_power_sum (:) = 0.0

  do k = lo(3), hi(3)
    do j = lo(2), hi(2)
      do i = lo(1), hi(1)

        ! Calculate grid points in k-space.
        if (i <= domain_grid_length/2) then
          k_index_i = i
        else
          k_index_i = i-domain_grid_length
        end if
        if (j <= domain_grid_length/2) then
          k_index_j = j
        else
          k_index_j = j-domain_grid_length
        end if
        if (k <= domain_grid_length/2) then
          k_index_k = k
        else
          k_index_k = k-domain_grid_length
        end if

        ! Skip 0 mode.
        if ((k_index_i == 0) .and. (k_index_j == 0) .and. (k_index_k == 0)) cycle

        ! Calculate magnitude of k vector.
        k_mag = k_fund * sqrt(real(k_index_i*k_index_i + k_index_j*k_index_j + k_index_k*k_index_k, c_double))

        ! Skip value if out of bounds.
        if ((k_mag < k_min) .or. (k_mag >= k_max)) cycle

        ! Find matching k bin.
        k_bin_index = -1
        do kk = 1, num_bins
        if (k_mag >= k_bin_edges(kk) .and. k_mag < k_bin_edges(kk+1)) then
            k_bin_index = kk
              exit
          end if
        end do
        if (k_bin_index < 0) stop 'COULD NOT MATCH K BIN'

        p_i = (overdensity_fft_real(i, j, k, 1)**2 + overdensity_fft_imag(i, j, k, 1)**2) / domain_length**3

        k_bin_count                (k_bin_index) = k_bin_count                (k_bin_index) + 1
        k_bin_power_weighted_k_sum (k_bin_index) = k_bin_power_weighted_k_sum (k_bin_index) + k_mag*p_i
        k_bin_power_sum            (k_bin_index) = k_bin_power_sum            (k_bin_index) + p_i

      end do
    end do
  end do

end subroutine calc_ps3d
