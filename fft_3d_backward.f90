subroutine fft_3d_backward (mf_fft_in_real, mf_fft_in_imag, lo, hi, domain_size, dx, comm, alloc_local, &
                            mf_fft_out_real, mf_fft_out_imag)
    use, intrinsic :: iso_c_binding
    use fftw3_mpi

    implicit none

    integer(c_int), intent(in ) :: lo(3), hi(3), domain_size(3), comm
    integer(c_intptr_t), intent(in) :: alloc_local
    real(c_double), intent(in) :: mf_fft_in_real (lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), 1)
    real(c_double), intent(in) :: mf_fft_in_imag (lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), 1)
    real(c_double), intent(in) :: dx

    real(c_double), intent(out) :: mf_fft_out_real(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), 1)
    real(c_double), intent(out) :: mf_fft_out_imag(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), 1)

    integer(c_intptr_t) :: l, m, n, local_n
    type(c_ptr) :: plan, cdata
    complex(c_double_complex), pointer :: mf_data(:, :, :)

    integer(c_int) :: cmplx_i, cmplx_j, cmplx_k
    integer(c_int) :: i, j, k
    real(c_double) :: grid_volume

    ! Since every process should have exactly 1 box, we won't need to
    ! make repeated calls to these init() functions.
    call fftw_mpi_init()

    l = domain_size(1)
    m = domain_size(2)
    n = domain_size(3)

    local_n = hi(3)-lo(3)+1

    ! We must use these special C-style pointers to feed data into FFTW.
    cdata = fftw_alloc_complex(alloc_local)
    call c_f_pointer(cdata, mf_data, [l, m, local_n])

    ! The sign convention used in the Nyx algorithms is the opposite of that used in FFTW. So a "forward" DFT in Nyx is "backward"
    ! in FFTW, and vice versa.
    plan = fftw_mpi_plan_dft_3d (n, m, l, mf_data, mf_data, comm, FFTW_FORWARD, FFTW_MEASURE)

    ! First translate the separate real and imaginary MultiFab data to a complex array. The indices for the local complex data
    ! arrays start their indices at one, unlike the Box indices which follow the domain grid coordinates.
    do k = lo(3), hi(3)
      cmplx_k = k - lo(3) + 1
      do j = lo(2), hi(2)
        cmplx_j = j - lo(2) + 1
        do i = lo(1), hi(1)
          cmplx_i = i - lo(1) + 1
          mf_data(cmplx_i, cmplx_j, cmplx_k) = cmplx(mf_fft_in_real(i, j, k, 1), mf_fft_in_imag(i, j, k, 1), c_double_complex)
        end do
      end do
    end do

    call fftw_mpi_execute_dft(plan, mf_data, mf_data)

    ! Copy data back to separate MultiFabs for real and imaginary components.
    do k = lo(3), hi(3)
      cmplx_k = k - lo(3) + 1
      do j = lo(2), hi(2)
        cmplx_j = j - lo(2) + 1
        do i = lo(1), hi(1)
          cmplx_i = i - lo(1) + 1
          ! Strip out real and imaginary components and save as separate
          ! MultiFabs. This syntax is super dumb.
          mf_fft_out_real (i, j, k, 1) = real(real (mf_data(cmplx_i, cmplx_j, cmplx_k)), c_double)
          mf_fft_out_imag (i, j, k, 1) = real(aimag(mf_data(cmplx_i, cmplx_j, cmplx_k)), c_double)
        end do
      end do
    end do

    call fftw_destroy_plan(plan)
    call fftw_free(cdata)

    ! Assumes cubic domain
    grid_volume = (domain_size(1) * dx)**3

    ! Normalize by grid volume
    mf_fft_out_real (:, :, :, 1) = mf_fft_out_real (:, :, :, 1) / grid_volume
    mf_fft_out_imag (:, :, :, 1) = mf_fft_out_imag (:, :, :, 1) / grid_volume

    call fftw_mpi_cleanup()

end subroutine fft_3d_backward
