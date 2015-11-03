! FFTW docs recommend doing this. They can't supply .mod files themselves because most .mod files are incompatible with different
! compilers, different compiler versions, etc.

module fftw3_mpi
    use, intrinsic :: iso_c_binding
    include 'fftw3-mpi.f03'
end module fftw3_mpi
