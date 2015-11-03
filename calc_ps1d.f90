! Calculate power spectrum P(k) of a 1-D pencil array of data. The "overdensity" input array is any scalar field (not just density)
! in the form g(x) = (x/x_bar) - 1, where x_bar is the mean.

subroutine calc_ps1d (overdensity_fft_real, overdensity_fft_imag, lo, hi, ng, dir, num_bins, k_bin_edges, domain_length, &
                 k_bin_count, k_bin_power_weighted_k_sum, k_bin_power_sum)

  use, intrinsic :: iso_c_binding
  implicit none

  integer(c_int), intent(in) :: lo(3), hi(3), ng, num_bins, dir
  real(c_double), intent(in), target :: overdensity_fft_real (lo(1)-ng:hi(1)+ng, lo(2)-ng:hi(2)+ng, lo(3)-ng:hi(3)+ng, 1)
  real(c_double), intent(in), target :: overdensity_fft_imag (lo(1)-ng:hi(1)+ng, lo(2)-ng:hi(2)+ng, lo(3)-ng:hi(3)+ng, 1)
  real(c_double), intent(in) :: k_bin_edges (num_bins+1)
  real(c_double), intent(in) :: domain_length

  integer(c_int), intent(out) :: k_bin_count (num_bins)
  real(c_double), intent(out) :: k_bin_power_weighted_k_sum (num_bins)
  real(c_double), intent(out) :: k_bin_power_sum (num_bins)

  integer(c_int) :: i
  integer(c_int) :: k
  integer(c_int) :: k_bin_index
  integer(c_int) :: k_index

  real(c_double) :: k_fund
  real(c_double) :: k_mag
  real(c_double) :: k_max
  real(c_double) :: k_min
  real(c_double) :: k_nyquist
  real(c_double) :: p_i
  real(c_double) :: pi

  real(c_double), pointer :: overdensity_fft_real_pencil (:)
  real(c_double), pointer :: overdensity_fft_imag_pencil (:)

  complex(c_double_complex), allocatable :: overdensity_fft_pencil (:)

  pi = acos(-1.0)

  ! Align 1-D array pointers along pencils to simplify indexing later on.
  if (dir == 0) then
      overdensity_fft_real_pencil => overdensity_fft_real (lo(1):hi(1), lo(2), lo(3), 1)
      overdensity_fft_imag_pencil => overdensity_fft_imag (lo(1):hi(1), lo(2), lo(3), 1)
  else if (dir == 1) then
      overdensity_fft_real_pencil => overdensity_fft_real (lo(1), lo(2):hi(2), lo(3), 1)
      overdensity_fft_imag_pencil => overdensity_fft_imag (lo(1), lo(2):hi(2), lo(3), 1)
  else if (dir == 2) then
      overdensity_fft_real_pencil => overdensity_fft_real (lo(1), lo(2), lo(3):hi(3), 1)
      overdensity_fft_imag_pencil => overdensity_fft_imag (lo(1), lo(2), lo(3):hi(3), 1)
  end if

  allocate(overdensity_fft_pencil(size(overdensity_fft_real_pencil)))

  ! Combine the real and imaginary parts into a single complex array.
  do i = 1, size(overdensity_fft_pencil)
    overdensity_fft_pencil(i) = cmplx(overdensity_fft_real_pencil(i), overdensity_fft_imag_pencil(i), c_double_complex)
  end do

  k_fund = 2.0 * pi / domain_length

  k_nyquist = k_fund * size(overdensity_fft_pencil)/2

  k_min = minval(k_bin_edges)
  k_max = maxval(k_bin_edges)

  if (k_min < 0.0 .or. k_max > k_nyquist) then
    print *, "(k_min, k_max, k_nyquist) = ", k_min, k_max, k_nyquist
    stop 'spectral k values out of bounds'
  end if

  k_bin_count (:) = 0
  k_bin_power_weighted_k_sum (:) = 0.0
  k_bin_power_sum (:) = 0.0

  do i = lbound(overdensity_fft_pencil, 1), ubound(overdensity_fft_pencil, 1)
    if (i <= size(overdensity_fft_pencil)/2+1) then
      k_index = i-1
    else
      k_index = i-1 - size(overdensity_fft_pencil)
    end if

    if (k_index == 0) cycle

    k_mag = k_fund * real(abs(k_index), c_double)

    if (k_mag < k_min .or. k_mag >= k_max) cycle

    ! Find the k bin.
    k_bin_index = -1
    do k = 1, num_bins
    if (k_mag >= k_bin_edges(k) .and. k_mag < k_bin_edges(k+1)) then
        k_bin_index = k
          exit
      end if
    end do
    if (k_bin_index < 0) stop 'COULD NOT MATCH K BIN'

    p_i = (abs(overdensity_fft_pencil(i)))**2 / domain_length

    k_bin_count                (k_bin_index) = k_bin_count                (k_bin_index) + 1
    k_bin_power_weighted_k_sum (k_bin_index) = k_bin_power_weighted_k_sum (k_bin_index) + k_mag*p_i
    k_bin_power_sum            (k_bin_index) = k_bin_power_sum            (k_bin_index) + p_i

  end do

  deallocate(overdensity_fft_pencil)

end subroutine calc_ps1d
