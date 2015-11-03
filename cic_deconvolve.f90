subroutine cic_deconvolve(dm_density_real, dm_density_imag, lo, hi, domain_size)
  use, intrinsic :: iso_c_binding
  implicit none

  integer(c_int), intent(in) :: lo(3), hi(3), domain_size(3)
  real(c_double), intent(inout) :: dm_density_real (lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), 1)
  real(c_double), intent(inout) :: dm_density_imag (lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), 1)

  integer(c_int) :: ki_nyquist_x, ki_nyquist_y, ki_nyquist_z

  real(c_double) :: pi
  real(c_double) :: sinc_arg_x, sinc_arg_y, sinc_arg_z
  real(c_double) :: sinc_arg_factor_x, sinc_arg_factor_y, sinc_arg_factor_z
  real(c_double) :: sinc_x, sinc_y, sinc_z

  integer(c_int) :: i, j, k
  integer(c_int) :: ki_x, ki_y, ki_z
  real(c_double) :: sinc_prod, cic_ft

  pi = acos(-1.0)

  ! Copied straight from old gimlet. Will this still work if the grid
  ! size is odd?
  ki_nyquist_x = domain_size(1) / 2
  ki_nyquist_y = domain_size(2) / 2
  ki_nyquist_z = domain_size(3) / 2

  sinc_arg_factor_x = pi / (2.0 * real(ki_nyquist_x, c_double))
  sinc_arg_factor_y = pi / (2.0 * real(ki_nyquist_y, c_double))
  sinc_arg_factor_z = pi / (2.0 * real(ki_nyquist_z, c_double))

  do k = lo(3), hi(3)

    if (k <= ki_nyquist_z) then
      ki_z = k
    else
      ki_z = domain_size(3) - k
    end if
    if (ki_z == 0) then
      sinc_z = 1.0
    else
      sinc_arg_z = ki_z * sinc_arg_factor_z
      sinc_z = sin(sinc_arg_z) / sinc_arg_z
    end if

    do j = lo(2), hi(2)

      if (j <= ki_nyquist_y) then
        ki_y = j
      else
        ki_y = domain_size(2) - j
      end if
      if (ki_y == 0) then
        sinc_y = 1.0
      else
        sinc_arg_y = ki_y * sinc_arg_factor_y
        sinc_y = sin(sinc_arg_y) / sinc_arg_y
      end if

      do i = lo(1), hi(1)

        if (i <= ki_nyquist_x) then
          ki_x = i
        else
          ki_x = domain_size(1) - i
        end if
        if (ki_x == 0) then
          sinc_x = 1.0
        else
          sinc_arg_x = ki_x * sinc_arg_factor_x
          sinc_x = sin(sinc_arg_x) / sinc_arg_x
        end if

        sinc_prod = sinc_x * sinc_y * sinc_z
        cic_ft = sinc_prod**2

        dm_density_real (i, j, k, 1) = dm_density_real (i, j, k, 1) / cic_ft
        dm_density_imag (i, j, k, 1) = dm_density_imag (i, j, k, 1) / cic_ft
      end do
    end do
  end do

end subroutine cic_deconvolve
