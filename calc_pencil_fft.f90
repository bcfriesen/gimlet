! Calculate 1-D Fourier transform of grid data along each pencil spanning the grid using FFTW. We can always do serial FFTs because
! the memory footprint for pencils is small even for very large grids.

subroutine calc_pencil_fft(overdensity, lo, hi, ng, dir, domain_length, dx, num_bins, k_bin_edges, &
                           overdensity_fft_real, overdensity_fft_imag)
    use, intrinsic :: iso_c_binding
    use fftw3
    implicit none

    integer(c_int), intent(in) :: lo(3), hi(3), ng
    real(c_double), intent(in), target :: overdensity (lo(1)-ng:hi(1)+ng, lo(2)-ng:hi(2)+ng, lo(3)-ng:hi(3)+ng, 1)
    integer(c_int), intent(in) :: dir
    real(c_double), intent(in) :: domain_length, dx
    integer(c_int), intent(in) :: num_bins
    real(c_double), intent(in) :: k_bin_edges (num_bins+1)

    real(c_double), intent(out), target :: overdensity_fft_real (lo(1)-ng:hi(1)+ng, lo(2)-ng:hi(2)+ng, lo(3)-ng:hi(3)+ng, 1)
    real(c_double), intent(out), target :: overdensity_fft_imag (lo(1)-ng:hi(1)+ng, lo(2)-ng:hi(2)+ng, lo(3)-ng:hi(3)+ng, 1)

    integer(c_int) :: i

    type(c_ptr) :: plan

    real(c_double), pointer :: overdensity_data(:)
    complex(c_double_complex), allocatable :: overdensity_fft_pencil (:)
    complex(c_double_complex), allocatable :: overdensity_pencil(:)

    real(c_double), pointer :: overdensity_fft_pencil_real(:)
    real(c_double), pointer :: overdensity_fft_pencil_imag(:)

    ! Align 1-D array pointers along pencils to simplify indexing later in the calculation.
    if (dir == 0) then
        overdensity_data            => overdensity          (lo(1):hi(1), lo(2), lo(3), 1)
        overdensity_fft_pencil_real => overdensity_fft_real (lo(1):hi(1), lo(2), lo(3), 1)
        overdensity_fft_pencil_imag => overdensity_fft_imag (lo(1):hi(1), lo(2), lo(3), 1)
    else if (dir == 1) then
        overdensity_data            => overdensity          (lo(1), lo(2):hi(2), lo(3), 1)
        overdensity_fft_pencil_real => overdensity_fft_real (lo(1), lo(2):hi(2), lo(3), 1)
        overdensity_fft_pencil_imag => overdensity_fft_imag (lo(1), lo(2):hi(2), lo(3), 1)
    else if (dir == 2) then
        overdensity_data            => overdensity          (lo(1), lo(2), lo(3):hi(3), 1)
        overdensity_fft_pencil_real => overdensity_fft_real (lo(1), lo(2), lo(3):hi(3), 1)
        overdensity_fft_pencil_imag => overdensity_fft_imag (lo(1), lo(2), lo(3):hi(3), 1)
    end if

    allocate(overdensity_pencil(size(overdensity_data)))
    allocate(overdensity_fft_pencil(size(overdensity_data)))

    ! The input data is real so the imaginary component of the complex pencil array is zero.
    do i = 1, size(overdensity_pencil)
      overdensity_pencil(i) = cmplx(overdensity_data(i), 0.0, c_double_complex)
    end do

    plan = fftw_plan_dft_1d(size(overdensity_pencil), overdensity_pencil, overdensity_fft_pencil, FFTW_BACKWARD, FFTW_ESTIMATE)

    call fftw_execute_dft(plan, overdensity_pencil, overdensity_fft_pencil)

    call fftw_destroy_plan(plan)

    ! Normalization
    overdensity_fft_pencil (:) = overdensity_fft_pencil (:) * dx

    ! Strip out real and imaginary components and save as separate MultiFabs. This syntax is weird because of the dual meanings of
    ! "real".
    do i = 1, size(overdensity_fft_pencil)
      overdensity_fft_pencil_real (i) = real(real (overdensity_fft_pencil(i)), c_double)
      overdensity_fft_pencil_imag (i) = real(aimag(overdensity_fft_pencil(i)), c_double)
    end do

    deallocate(overdensity_pencil)
    deallocate(overdensity_fft_pencil)

end subroutine calc_pencil_fft
