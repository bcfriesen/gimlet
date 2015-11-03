! Calculate optical depth of a scalar field divided into pencils pointing along a particular axis.

subroutine calc_tau (state_data, eos_data, mom, &
                     lo, hi, ng_mom, ng_state, ng_eos, ng_tau, mean_density, dx, dir, z, domain_length, omega_m, omega_l, omega_b, &
                     h_0, tau)

    use, intrinsic :: iso_c_binding
    use meth_params_module, only: NVAR, UEINT, URHO, TEMP_COMP
    use eos_module, only: nyx_eos_nh0_and_nhep
    use fundamental_constants_module, only: e_to_cgs, density_to_cgs
    implicit none

    integer(c_int), intent(in) :: lo(3), hi(3), ng_mom, ng_state, ng_eos, ng_tau, dir
    real(c_double), intent(in), target :: state_data (lo(1)-ng_state:hi(1)+ng_state, &
                                                      lo(2)-ng_state:hi(2)+ng_state, &
                                                      lo(3)-ng_state:hi(3)+ng_state, NVAR)
    ! I guess we always hard-code the # of EOS components as 2?
    real(c_double), intent(in), target :: eos_data   (lo(1)-ng_eos:hi(1)+ng_eos, &
                                                      lo(2)-ng_eos:hi(2)+ng_eos, &
                                                      lo(3)-ng_eos:hi(3)+ng_eos, 2)
    real(c_double), intent(in), target :: mom        (lo(1)-ng_mom:hi(1)+ng_mom, &
                                                      lo(2)-ng_mom:hi(2)+ng_mom, &
                                                      lo(3)-ng_mom:hi(3)+ng_mom, 1)
    real(c_double), intent(in) :: mean_density, dx, z, domain_length, omega_m, omega_l, omega_b, h_0

    real(c_double), intent(out), target :: tau (lo(1)-ng_tau:hi(1)+ng_tau, &
                                                lo(2)-ng_tau:hi(2)+ng_tau, &
                                                lo(3)-ng_tau:hi(3)+ng_tau, 1)

    real(c_double), pointer :: tau_pencil (:)

    ! Pointers to pencils in specified directions through the above MultiFabs.
    real(c_double), pointer :: density_pencil (:)
    real(c_double), pointer :: e_int_pencil   (:)
    real(c_double), pointer :: temp_pencil    (:)
    real(c_double), pointer :: mom_pencil     (:)

    real(c_double), parameter :: k_cgs = 1.38064d-16
    real(c_double), parameter :: m_p_cgs = 1.672623d-24
    real(c_double), parameter :: hilya_sigma0_cgs = 0.011045900074845754d0
    real(c_double), parameter :: hilya_lambda_cgs = 1.2150227340678867d-05
    real(c_double), parameter :: mpc_cgs = 3.08568d24
    real(c_double), parameter :: v_thermal_factor = sqrt(2.0d0 * k_cgs / m_p_cgs)

    real(c_double) :: hubble_z
    real(c_double) :: tau_factor

    real(c_double) :: cell_nhi
    real(c_double) :: cell_vlos
    real(c_double) :: cell_rho

    real(c_double) :: cell_v_lo
    real(c_double) :: cell_v_mid
    real(c_double) :: cell_v_hi

    real(c_double) :: cell_v_size
    real(c_double) :: cell_vth
    real(c_double) :: cell_temp

    real(c_double) :: erf_arg_lo
    real(c_double) :: erf_arg_hi

    real(c_double) :: domain_velocity_size

    integer(c_int) :: i
    integer(c_int) :: ipix

    real(c_double) :: pixel_v_size
    real(c_double) :: tau_common
    real(c_double) :: v_doppler
    real(c_double) :: v_pixel

    real(c_double) :: v_lc
    real(c_double) :: pix_broad
    real(c_double) :: pix_lc

    real(c_double) :: comoving_a
    real(c_double) :: a3

    real(c_double) :: h

    real(c_double) :: cell_density_cgs
    real(c_double) :: rho_b_cgs
    real(c_double) :: e_int
    real(c_double) :: cell_nhep

    integer(c_int) :: ipix_lo
    integer(c_int) :: ipix_hi
    integer(c_int) :: ipix_wrapped
    integer(c_int) :: num_pixels
    integer(c_int) :: num_sig_pixels

    real(c_double), parameter :: rho_crit_100_cgs = 1.8788200387d-29
    real(c_double) :: rho_eos_units

    comoving_a = 1.0 / (1.0 + z)
    a3 = comoving_a * comoving_a * comoving_a

    ! Hubble parameter
    h = h_0 / 100.0

    hubble_z = h_0 * sqrt(omega_m * (1.0d0 + z)**3 + omega_l) * 1.0d5 / mpc_cgs
    tau_factor = hilya_sigma0_cgs * hilya_lambda_cgs / (2.0d0 * hubble_z)

    domain_velocity_size = 1.0d7 * domain_length * sqrt(omega_m * (1.0d0 + z)**3 + omega_l) / (1.0d0 + z)

    ! If we've called calc_tau() then it is assumed that we have already regridded all of the input MultiFabs into single-cell-thick
    ! pencils. (If this is not true, then the behavior of this routine will be undefined.) This means lo() = hi() for the 2
    ! dimensions perpendicular to the pencils. Below I have arbitrarily chosen lo() as the argument for the non-pencil dimensions,
    ! but could equivalently have chosen hi().

    ! Pencils point in x-direction.
    if (dir == 0) then
      density_pencil => state_data (lo(1):hi(1), lo(2), lo(3), URHO     )
      e_int_pencil   => state_data (lo(1):hi(1), lo(2), lo(3), UEINT    )
      temp_pencil    => eos_data   (lo(1):hi(1), lo(2), lo(3), TEMP_COMP)
      mom_pencil     => mom        (lo(1):hi(1), lo(2), lo(3), 1        )
      tau_pencil     => tau        (lo(1):hi(1), lo(2), lo(3), 1        )
    ! Pencils point in y-direction.
    else if (dir == 1) then
      density_pencil => state_data (lo(1), lo(2):hi(2), lo(3), URHO     )
      e_int_pencil   => state_data (lo(1), lo(2):hi(2), lo(3), UEINT    )
      temp_pencil    => eos_data   (lo(1), lo(2):hi(2), lo(3), TEMP_COMP)
      mom_pencil     => mom        (lo(1), lo(2):hi(2), lo(3), 1        )
      tau_pencil     => tau        (lo(1), lo(2):hi(2), lo(3), 1        )
    ! Pencils point in z-direction.
    else if (dir == 2) then
      density_pencil => state_data (lo(1), lo(2), lo(3):hi(3), URHO     )
      e_int_pencil   => state_data (lo(1), lo(2), lo(3):hi(3), UEINT    )
      temp_pencil    => eos_data   (lo(1), lo(2), lo(3):hi(3), TEMP_COMP)
      mom_pencil     => mom        (lo(1), lo(2), lo(3):hi(3), 1        )
      tau_pencil     => tau        (lo(1), lo(2), lo(3):hi(3), 1        )
    end if

    cell_v_size = domain_velocity_size / real(size(tau_pencil), c_double)
    pixel_v_size = domain_velocity_size / real(size(tau_pencil), c_double)

    tau_pencil(:) = 0.0

    num_pixels = size(tau_pencil)

    ! Behold a confusing feature in Fortran: the indices of an array pointer always use standard indexing, i.e., start at 1, even if
    ! it points to a target array with non-standard indexing. So I can't use "do i = lo(), hi()" when looping over these array
    ! pointers even though that's the typical approach in BoxLib Fortran functions.

    do i = lbound(tau_pencil, 1), ubound(tau_pencil, 1)

      cell_rho = density_pencil(i)
      e_int = e_int_pencil(i) / cell_rho

      rho_eos_units = cell_rho * omega_b * h*h * rho_crit_100_cgs / (a3 * mean_density)

      ! Get H and He densities from the EOS.
      call nyx_eos_nh0_and_nhep(z, rho_eos_units, e_int*e_to_cgs, cell_nhi, cell_nhep)

      cell_vlos = mom_pencil(i) / density_pencil(i) * 1.0e5

      cell_temp = temp_pencil(i)

      ! Thermal velocity
      cell_vth = v_thermal_factor * sqrt(cell_temp)
      v_doppler = cell_vth

      ! Cell velocities
      cell_v_lo  = cell_v_size * (real(i, c_double) - 1.0)
      cell_v_mid = cell_v_size * (real(i, c_double) - 0.5)
      cell_v_hi  = cell_v_size * (real(i, c_double)      )

      tau_common = tau_factor * cell_nhi

      v_lc = cell_v_mid + cell_vlos

      pix_lc = v_lc / pixel_v_size
      pix_broad = v_doppler / pixel_v_size

      num_sig_pixels = int(4.0 * pix_broad + 0.5, c_int) + 1

      ipix_lo = pix_lc - num_sig_pixels - 1
      ipix_hi = pix_lc + num_sig_pixels

      do ipix = ipix_lo, ipix_hi

        v_pixel = pixel_v_size * (real(ipix, c_double) + 0.5)

        erf_arg_hi = (v_pixel - cell_vlos - cell_v_hi) / v_doppler
        erf_arg_lo = (v_pixel - cell_vlos - cell_v_lo) / v_doppler

        ! Accounts for periodic boundary conditions. If the contributed pixel in index space is off the end of the skewer, wrap it
        ! around on the other end.
        ipix_wrapped = ipix - num_pixels * int(real(ipix, c_double)/real(num_pixels, c_double), c_int)
        if (ipix_wrapped < 1) then
          ipix_wrapped = ipix_wrapped + num_pixels
        end if

        tau_pencil(ipix_wrapped) = tau_pencil(ipix_wrapped) + tau_common * (erf(erf_arg_lo) - erf(erf_arg_hi))
      end do

    end do

end subroutine calc_tau
