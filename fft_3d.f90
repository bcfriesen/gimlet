! Calculate 3-D Fourier transform of a scalar field using distributed
! memory (MPI) version of FFTW. This routine requires that the Box
! distribution of the input MultiFab already conforms to that required
! by FFTW, i.e., striped along last (z) dimension.

subroutine fft_3d(mf_fft_in, lo, hi, domain_size, dx, comm, alloc_local, mf_fft_out_real, mf_fft_out_imag, threads_ok)
    use, intrinsic :: iso_c_binding
    use fftw3_mpi
    use omp_lib

    implicit none

    integer(c_int), intent(in ) :: lo(3), hi(3), domain_size(3), comm
    integer(c_intptr_t), intent(in) :: alloc_local
    real(c_double), intent(in) :: mf_fft_in (lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), 1)
    real(c_double), intent(in) :: dx
    integer, intent(in) :: threads_ok

    real(c_double), intent(out) :: mf_fft_out_real(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), 1)
    real(c_double), intent(out) :: mf_fft_out_imag(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), 1)

    integer(c_intptr_t) :: l, m, n, local_n
    type(c_ptr) :: plan, cdata
    complex(c_double_complex), pointer :: mf_data(:, :, :)

    integer :: cmplx_i, cmplx_j, cmplx_k
    integer :: i, j, k
    real(c_double) :: dx3

    dx3 = dx*dx*dx

    l = domain_size(1)
    m = domain_size(2)
    n = domain_size(3)

    local_n = hi(3)-lo(3)+1

    ! We must use these special C-style pointers to feed data into FFTW.
    cdata = fftw_alloc_complex(alloc_local)
    call c_f_pointer(cdata, mf_data, [l, m, local_n])

    if (threads_ok /= 0) call fftw_plan_with_nthreads(omp_get_num_threads())
    plan = fftw_mpi_plan_dft_3d (n, m, l, mf_data, mf_data, comm, FFTW_BACKWARD, FFTW_MEASURE)

    ! First translate the real MultiFab data to a complex array. The
    ! indices for the local complex data arrays start their indices at
    ! one, unlike the Box indices which follow the domain grid
    ! coordinates.
    do k = lo(3), hi(3)
      cmplx_k = k - lo(3) + 1
      do j = lo(2), hi(2)
        cmplx_j = j - lo(2) + 1
        do i = lo(1), hi(1)
          cmplx_i = i - lo(1) + 1
          ! Input data to the DFT is real data so the imaginary
          ! component is zero.
          mf_data(cmplx_i, cmplx_j, cmplx_k) = cmplx(mf_fft_in(i, j, k, 1), 0.0, c_double_complex)
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
          ! MultiFabs. This syntax is weird due to the dual meanings of
          ! "real".
          mf_fft_out_real (i, j, k, 1) = real(real (mf_data(cmplx_i, cmplx_j, cmplx_k)), c_double)
          mf_fft_out_imag (i, j, k, 1) = real(aimag(mf_data(cmplx_i, cmplx_j, cmplx_k)), c_double)
        end do
      end do
    end do

    call fftw_destroy_plan(plan)
    call fftw_free(cdata)

    ! Normalize by cell volume
    mf_fft_out_real (:, :, :, 1) = mf_fft_out_real (:, :, :, 1) * dx3
    mf_fft_out_imag (:, :, :, 1) = mf_fft_out_imag (:, :, :, 1) * dx3

end subroutine fft_3d
