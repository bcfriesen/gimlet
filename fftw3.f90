! FFTW docs recommend doing this. They can't supply .mod files themselves because most .mod files are incompatible with different
! compilers, different compiler versions, etc.

module fftw3
    use, intrinsic :: iso_c_binding
    include 'fftw3.f03'
end module fftw3
