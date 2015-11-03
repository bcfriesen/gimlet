! Calculate P(k, mu) power spectrum from Fourier transformed data along pencils. The input data here is called "overdensity" but
! it's actually any scalar field in the form f(x) = (x/x_bar) - 1 where x_bar is the mean value.

subroutine calc_k_mu_ps1d (overdensity_fft_real, overdensity_fft_imag, &
                           lo, hi, ng, dir, domain_length, &
                           num_k_bins, k_bin_edges, &
                           num_mu_bins, mu_bin_edges, &
                           bin_count, power_weighted_k_sum, power_weighted_mu_sum, power_sum)

  use, intrinsic :: iso_c_binding
  implicit none

  integer(c_int), intent(in) :: lo(3), hi(3), ng, num_k_bins, num_mu_bins, dir
  real(c_double), intent(in), target :: overdensity_fft_real (lo(1)-ng:hi(1)+ng, lo(2)-ng:hi(2)+ng, lo(3)-ng:hi(3)+ng, 1)
  real(c_double), intent(in), target :: overdensity_fft_imag (lo(1)-ng:hi(1)+ng, lo(2)-ng:hi(2)+ng, lo(3)-ng:hi(3)+ng, 1)
  real(c_double), intent(in) :: k_bin_edges (num_k_bins+1)
  real(c_double), intent(in) :: mu_bin_edges (num_mu_bins+1)
  real(c_double), intent(in) :: domain_length

  integer(c_int), intent(out) :: bin_count (num_k_bins*num_mu_bins)
  real(c_double), intent(out) :: power_weighted_k_sum (num_k_bins*num_mu_bins)
  real(c_double), intent(out) :: power_weighted_mu_sum (num_k_bins*num_mu_bins)
  real(c_double), intent(out) :: power_sum (num_k_bins*num_mu_bins)

  integer(c_int) :: i
  integer(c_int) :: k
  integer(c_int) :: k_bin_index
  integer(c_int) :: mu_bin_index
  integer(c_int) :: i_bin
  integer(c_int) :: k_index(3)

  real(c_double) :: k_fund
  real(c_double) :: k_mag
  real(c_double) :: k_max
  real(c_double) :: mu_min
  real(c_double) :: mu_max
  real(c_double) :: k_min
  real(c_double) :: k_nyquist
  real(c_double) :: p_i
  real(c_double) :: pi
  real(c_double) :: k_para
  real(c_double) :: mu

  real(c_double), pointer :: overdensity_fft_real_pencil (:)
  real(c_double), pointer :: overdensity_fft_imag_pencil (:)

  complex(c_double_complex), allocatable :: overdensity_fft_pencil (:)

  pi = acos(-1.0)

  ! Align array pointers along pencil direction. We use 1-D pointers here in order to simplify the indices in the calculations
  ! below.
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

  mu_min = minval(mu_bin_edges)
  mu_max = maxval(mu_bin_edges)

  if (k_min < 0.0 .or. k_max > k_nyquist) then
    print *, "(k_min, k_max, k_nyquist) = ", k_min, k_max, k_nyquist
    stop
  end if

  bin_count (:) = 0
  power_weighted_k_sum (:) = 0.0
  power_weighted_mu_sum (:) = 0.0
  power_sum (:) = 0.0

  do i = lbound(overdensity_fft_pencil, 1), ubound(overdensity_fft_pencil, 1)
    if (i <= size(overdensity_fft_pencil)/2+1) then
      k_index(dir+1) = i-1
    else
      k_index(dir+1) = i-1 - size(overdensity_fft_pencil)
    end if

    ! Assign wavenumber indices from the non-pencil grid directions. For the 2 non-pencil-pointing directions, lo() = hi() (since
    ! the pencils are only 1 cell thick) so there's nothing special about choosing lo() here.
    if (dir == 0) then
      if (lo(2) <= size(overdensity_fft_pencil)/2) then
        k_index(2) = lo(2)
      else
        k_index(2) = lo(2) - size(overdensity_fft_pencil)
      end if
      if (lo(3) <= size(overdensity_fft_pencil)/2) then
        k_index(3) = lo(3)
      else
        k_index(3) = lo(3) - size(overdensity_fft_pencil)
      end if
    else if (dir == 1) then
      if (lo(1) <= size(overdensity_fft_pencil)/2) then
        k_index(1) = lo(1)
      else
        k_index(1) = lo(1) - size(overdensity_fft_pencil)
      end if
      if (lo(3) <= size(overdensity_fft_pencil)/2) then
        k_index(3) = lo(3)
      else
        k_index(3) = lo(3) - size(overdensity_fft_pencil)
      end if
    else if (dir == 2) then
      if (lo(1) <= size(overdensity_fft_pencil)/2) then
        k_index(1) = lo(1)
      else
        k_index(1) = lo(1) - size(overdensity_fft_pencil)
      end if
      if (lo(2) <= size(overdensity_fft_pencil)/2) then
        k_index(2) = lo(2)
      else
        k_index(2) = lo(2) - size(overdensity_fft_pencil)
      end if
    end if

    ! Skip the zero mode.
    if (k_index(1) == 0 .and. k_index(2) == 0 .and. k_index(3) == 0) cycle

    k_para = k_fund * real(abs(k_index(dir+1)), c_double)
    ! Get magnitude of k vector.
    k_mag = k_fund * sqrt(real(k_index(1)*k_index(1) + k_index(2)*k_index(2) + k_index(3)*k_index(3), c_double))
    mu = k_para / k_mag

    ! Skip wavenumbers and angles outside the requested range.
    if (k_mag < k_min .or. k_mag >= k_max .or. mu < mu_min .or. mu > mu_max) cycle

    ! Match k bin.
    k_bin_index = -1
    do k = 1, num_k_bins
      if (k_mag >= k_bin_edges(k) .and. k_mag < k_bin_edges(k+1)) then
        k_bin_index = k
        exit
      end if
    end do
    if (k_bin_index < 0) stop 'COULD NOT MATCH K BIN'

    ! Match mu bin.
    mu_bin_index = -1
    do k = 1, num_mu_bins
      if (mu >= mu_bin_edges(k) .and. mu < mu_bin_edges(k+1)) then
        mu_bin_index = k
        exit
      end if
    end do
    if (mu == mu_max) mu_bin_index = num_mu_bins
    if (mu_bin_index < 0) stop 'COULD NOT MATCH MU BIN'

    ! Flatten index for (k, mu) arrays.
    i_bin = (k_bin_index-1)*num_mu_bins + mu_bin_index

    ! VOLUME DENOMINATOR ASSUMES CUBIC DOMAIN
    p_i = (abs(overdensity_fft_pencil(i)))**2 / domain_length**3

    bin_count             (i_bin) = bin_count             (i_bin) + 1
    power_weighted_k_sum  (i_bin) = power_weighted_k_sum  (i_bin) + p_i * k_mag
    power_weighted_mu_sum (i_bin) = power_weighted_mu_sum (i_bin) + p_i * mu
    power_sum             (i_bin) = power_sum             (i_bin) + p_i

  end do

  deallocate(overdensity_fft_pencil)

end subroutine calc_k_mu_ps1d
