! Query FFTW to get the Box distribution it prefers, given the total
! size of the grid. We'll use these results to construct a new BoxArray
! onto which to regrid our MultiFabs before we actually do the FFT.

subroutine get_fftw_box_sizes(domain_size, comm, local_z, local_k_offset, alloc_local)
  use, intrinsic :: iso_c_binding
  use fftw3_mpi
  implicit none

  integer(c_intptr_t), intent(in) :: domain_size(3)
  integer(c_int), intent(in) :: comm

  integer(c_intptr_t), intent(out) :: local_z, local_k_offset, alloc_local

  call fftw_mpi_init()

  ! Since we'll be doing the FFT in Fortran, we need to use
  ! Fortran-style array indexing, so we reverse the array size ordering
  ! to FFTW functions. (See FFTW docs for details.)
  alloc_local = fftw_mpi_local_size_3d (domain_size(3), domain_size(2), domain_size(1), comm, local_z, local_k_offset)

  call fftw_mpi_cleanup()

end subroutine get_fftw_box_sizes
