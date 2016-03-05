#include <string>
#include <new>
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <list>
#include <map>
#include <algorithm>
#include <fstream>
#include <sstream>

#include <MultiFab.H>
#include <ParallelDescriptor.H>
#include <Utility.H>
#include <FArrayBox.H>
#include <Geometry.H>
#include <BLFort.H>

#include <MakePencilBoxes.H>
#include <postprocess_tau_fields.H>
#include <pdf.H>
#include <ps3d.H>
#include <temperature_density_pdf2d.H>
#include <MakeFFTWBoxes.H>

BL_FORT_PROC_DECL(CALC_TAU, calc_tau) (
   const Real *density,
   const Real *temperature,
   const Real *e_int,
   const Real *mom,
   const int  *lo,
   const int  *hi,
   const int  *ng_mom,
   const int  *ng_state,
   const int  *ng_eos,
   const int  *ng_tau,
   const Real *mean_density,
   const Real *dx,
   const int  *dir,
   const Real *z,
   const Real *domain_length,
   const Real *omega_m,
   const Real *omega_l,
   const Real *omega_b,
   const Real *H_0,
   Real       *tau);

BL_FORT_PROC_DECL(CALC_RHO_M, calc_rho_m) (
  const Real* rho_b,
  const Real* rho_dm,
  const int* lo,
  const int* hi,
  const int* density_num_ghosts,
  const int* dm_density_num_ghosts,
  const int* rho_m_num_ghosts,
  const Real* omega_b,
  const Real* omega_m,
  const Real* mean_density,
  const Real* mean_dm_density,
  const Real* rho_m);

BL_FORT_PROC_DECL(CALC_N_HI, calc_n_hi) (
    const Real* z,
    const Real* density,
    const Real* e_int,
    const Real* mean_density,
    const Real* omega_b,
    const Real* H_0,
    const int* lo,
    const int* hi,
    const int* state_ng,
    const int* n_hi_ng,
    const Real* n_hi);

BL_FORT_PROC_DECL(CALC_ABS_V, calc_abs_v) (
    const Real* xmom,
    const Real* ymom,
    const Real* zmom,
    const int* mom_num_ghosts,
    const Real* density,
    const int* density_num_ghosts,
    const int* lo,
    const int* hi,
    const Real* abs_v);

BL_FORT_PROC_DECL(CALC_ABS_VZ, calc_abs_vz) (
    const Real* zmom,
    const Real* density,
    const int* zmom_num_ghosts,
    const int* density_num_ghosts,
    const int* vz_num_ghosts,
    const int* lo,
    const int* hi,
    const Real* abs_vz);

BL_FORT_PROC_DECL(FFT_3D, fft_3d) (
        const Real*     mf_fft_in,
        const int*      lo,
        const int*      hi,
        const int*      domain_size_int,
        const Real*     dx,
        const int*      comm,
        const intptr_t* alloc_local,
        const Real*     mf_fft_out_real,
        const Real*     mf_fft_out_imag);

BL_FORT_PROC_DECL(FFT_3D_BACKWARD, fft_3d_backward) (
        const Real*     mf_fft_in_real,
        const Real*     mf_fft_in_imag,
        const int*      lo,
        const int*      hi,
        const int*      domain_size_int,
        const Real*     dx,
        const int*      comm,
        const intptr_t* alloc_local,
        const Real*     mf_fft_out_real,
        const Real*     mf_fft_out_imag);

BL_FORT_PROC_DECL(CIC_DECONVOLVE, cic_deconvolve) (
    const Real* dm_density_real,
    const Real* dm_density_imag,
    const int* lo,
    const int* hi,
    const int* domain_size);


void
do_analysis(const Real     omega_b,
            const Real     omega_m,
            const Real     omega_l,
            const Real     h,
            const Real     comoving_a,
            const MultiFab &density,
            const MultiFab &temperature,
            const MultiFab &e_int,
            const MultiFab &dm_density,
            const MultiFab &xmom,
            const MultiFab &ymom,
            const MultiFab &zmom,
            const Geometry &geom,
            const int      nStep)
{

    const Real z = (1.0/comoving_a) - 1.0;

    if (ParallelDescriptor::IOProcessor()) {
        // Undo the effects of std::scientific, since these numbers are usually
        // of order unity.
        std::cout.setf(std::ios::fixed, std::ios::floatfield);
        std::cout.setf(std::ios::showpoint);
        std::cout << "===== gimlet got these cosmological parameters: " << std::endl;
        std::cout << "=====      Omega_B = " << omega_b    << std::endl;
        std::cout << "=====      Omega_M = " << omega_m    << std::endl;
        std::cout << "=====      Omega_L = " << omega_l    << std::endl;
        std::cout << "=====      h       = " << h          << std::endl;
        std::cout << "=====      a       = " << comoving_a << std::endl;
        std::cout << "=====      z       = " << z          << std::endl;
        std::cout << std::flush;
    }

    std::cout << std::scientific;

    const Box& problem_domain = geom.Domain();
    const long num_cells = problem_domain.numPts();

    const Real mean_density = density.norm1() / Real(num_cells);
    if (ParallelDescriptor::IOProcessor())
        std::cout << "===== Mean baryon density = " << std::scientific << mean_density << " Msun/Mpc^3" << std::endl;

    const BoxArray& ba1 = density.boxArray();

    /* Every box will be a one-cell-thick pencil. We will trace rays in all 3
     * dimensions, so we will re-grid 3 times. The boxes will be (long) x 1 x 1
     * or 1 x (long) x 1 or 1 x 1 x (long). */

    const Real domain_length = geom.ProbLength(0); // units are Mpc

    const Real H_0 = h*100.0;

    Real time1, total_time;

    for (int dir = 0; dir < BL_SPACEDIM; ++dir)
    {
        if (ParallelDescriptor::IOProcessor())
          std::cout << std::setfill('=') << std::setw(40) << " Copying data to pencils in dir = " << dir << " ... " << std::flush;

        time1 = ParallelDescriptor::second();
        BoxList box_list = MakePencilBoxes(geom, dir);
        BoxArray ba(box_list);

        MultiFab density_pencils     (ba,     density.nComp(),     density.nGrow());
        MultiFab temperature_pencils (ba, temperature.nComp(), temperature.nGrow());
        MultiFab e_int_pencils       (ba,       e_int.nComp(),       e_int.nGrow());

        int num_ghosts_mom;
        if (dir == 0) {
            num_ghosts_mom = xmom.nGrow();
        } else if (dir == 1) {
            num_ghosts_mom = ymom.nGrow();
        } else if (dir == 2){
            num_ghosts_mom = zmom.nGrow();
        }
        const int num_ghosts_state_data =     density.nGrow();
        const int num_ghosts_eos_data   = temperature.nGrow();

        MultiFab mom_pencils (ba, 1, num_ghosts_mom);

        density_pencils.copy(density);
        temperature_pencils.copy(temperature);
        e_int_pencils.copy(e_int);
        if (dir == 0) {
            mom_pencils.copy(xmom);
        } else if (dir == 1) {
            mom_pencils.copy(ymom);
        } else if (dir == 2) {
            mom_pencils.copy(zmom);
        }

        total_time = ParallelDescriptor::second() - time1;
        if (ParallelDescriptor::IOProcessor())
          std::cout << std::setw(15) << " done. (" << total_time << " sec)" << std::endl;

        const int num_ghosts_tau = 0;
        MultiFab tau_redshift_space_pencils (ba, 1, num_ghosts_tau);

        // mesh spacing
        Array<Real> dx(BL_SPACEDIM);
        for (unsigned int n = 0; n < BL_SPACEDIM; ++n)
          dx.at(n) = (geom.ProbHi(n) - geom.ProbLo(n))/problem_domain.length(n);

        if (ParallelDescriptor::IOProcessor())
          std::cout << std::setfill('=') << std::setw(46) << " Calculating optical depth ... " << std::flush;

        time1 = ParallelDescriptor::second();
        for (MFIter mfi(density_pencils); mfi.isValid(); ++mfi) {
          const Box& bx = mfi.validbox();
          BL_FORT_PROC_CALL(CALC_TAU, calc_tau)
              ((density_pencils)[mfi].dataPtr(),
              (temperature_pencils)[mfi].dataPtr(),
              (e_int_pencils)[mfi].dataPtr(),
              (mom_pencils)[mfi].dataPtr(),
              bx.loVect(),
              bx.hiVect(),
              &num_ghosts_mom,
              &num_ghosts_state_data,
              &num_ghosts_eos_data,
              &num_ghosts_tau,
              &mean_density,
              &(dx[0]),
              &dir,
              &z,
              &domain_length,
              &omega_m,
              &omega_l,
              &omega_b,
              &H_0,
              (tau_redshift_space_pencils)[mfi].dataPtr());
        }
        total_time = ParallelDescriptor::second() - time1;

        if (ParallelDescriptor::IOProcessor())
          std::cout << std::setw(15) << " done. (" << total_time << " sec)" << std::endl;

        // Copy pencil data back to normally shaped boxes so other
        // post-processing tasks (which don't depend on pencils) will go
        // faster.
        MultiFab tau_redshift_space(ba1, 1, num_ghosts_tau);
        tau_redshift_space.copy(tau_redshift_space_pencils);

        postprocess_tau_fields(tau_redshift_space, geom, dir, nStep);
    }

    // These post-processing steps don't depend on pencils or LOS directions,
    // so we only do them once.

    // Braces force these temporary variables to go out of scope as soon as
    // we're done with them. Otherwise they linger until the end of the
    // analysis routine, which can make its memory footprint very large.

    {
        if (ParallelDescriptor::IOProcessor())
            std::cout << std::endl << std::setfill('=') << std::setw(46) << " Loading rho_b data ... " << std::flush;
        time1 = ParallelDescriptor::second();
        // Density PDF and power spectrum need to be in mean units
        const Real mean_density = density.norm1() / Real(num_cells);
        MultiFab density_divided_by_mean(density.boxArray(), 1, 0);
        density_divided_by_mean.copy(density);
        density_divided_by_mean.mult(1.0/mean_density);
        total_time = ParallelDescriptor::second() - time1;
        if (ParallelDescriptor::IOProcessor())
          std::cout << std::setw(15) << " done. (" << total_time << " sec)" << std::endl;

        if (ParallelDescriptor::IOProcessor())
          std::cout << std::setfill('=') << std::setw(46) << " Calculating rho_b PDF ... " << std::flush;
        #warning TODO: make density PDF knobs run-time parameters
        time1 = ParallelDescriptor::second();
        pdf (density_divided_by_mean, geom, nStep, "rhob", -2.0, 5.0, 400);
        total_time = ParallelDescriptor::second() - time1;
        if (ParallelDescriptor::IOProcessor())
          std::cout << std::setw(15) << " done. (" << total_time << " sec)" << std::endl;

        if (ParallelDescriptor::IOProcessor())
            std::cout << std::setfill('=') << std::setw(46) << " Calculating rho_b P(k) ... " << std::flush;
        #warning TODO: make density power spectrum knobs run-time parameters
        time1 = ParallelDescriptor::second();
        ps3d (density, geom, nStep, "rhob");
        total_time = ParallelDescriptor::second() - time1;
        if (ParallelDescriptor::IOProcessor())
            std::cout << std::setfill('=') << std::setw(15) << " done. (" << total_time << " sec)" << std::endl;
    }

    {
        if (ParallelDescriptor::IOProcessor())
          std::cout << std::setfill('=') << std::setw(46) << " Calculating T PDF ... " << std::flush;
        #warning TODO: make temperature PDF knobs run-time parameters
        time1 = ParallelDescriptor::second();
        pdf (temperature, geom, nStep, "temp", 3.0, 8.0, 400);
        total_time = ParallelDescriptor::second() - time1;
        if (ParallelDescriptor::IOProcessor())
          std::cout << std::setw(15) << " done. (" << total_time << " sec)" << std::endl;

        if (ParallelDescriptor::IOProcessor())
          std::cout << std::setfill('=') << std::setw(46) << " Calculating T P(k) ... " << std::flush;
        #warning TODO: make temperature power spectrum knobs run-time parameters
        time1 = ParallelDescriptor::second();
        ps3d (temperature, geom, nStep, "temp");
        total_time = ParallelDescriptor::second() - time1;
        if (ParallelDescriptor::IOProcessor())
            std::cout << std::setw(15) << " done. (" << total_time << " sec)" << std::endl;
    }

    {
        const Real mean_density = density.norm1() / Real(num_cells);
        if (ParallelDescriptor::IOProcessor())
          std::cout << std::setfill('=') << std::setw(46) << " Calculating T-rho_b 2-D PDF ... " << std::flush;
        #warning TODO: make temperature PDF knobs run-time parameters
        time1 = ParallelDescriptor::second();
        temperature_density_pdf2d (temperature, density, geom, mean_density, nStep);
        total_time = ParallelDescriptor::second() - time1;
        if (ParallelDescriptor::IOProcessor())
            std::cout << std::setw(15) << " done. (" << total_time << " sec)" << std::endl;
    }

    {
        if (ParallelDescriptor::IOProcessor())
            std::cout << std::endl << std::setfill('=') << std::setw(46) << " Loading rho_dm data ... " << std::flush;
        time1 = ParallelDescriptor::second();
        const Real mean_dm_density = dm_density.norm1() / Real(num_cells);
        MultiFab dm_density_divided_by_mean(ba1, 1, 0);
        dm_density_divided_by_mean.copy(dm_density);
        dm_density_divided_by_mean.mult(1.0/mean_dm_density);
        total_time = ParallelDescriptor::second() - time1;
        if (ParallelDescriptor::IOProcessor())
            std::cout << std::setw(15) << " done. (" << total_time << " sec)" << std::endl;

        if (ParallelDescriptor::IOProcessor())
          std::cout << std::setfill('=') << std::setw(46) << " Calculating rho_dm PDF ... " << std::flush;
        time1 = ParallelDescriptor::second();
        #warning TODO: make dm_density PDF knobs run-time parameters
        pdf (dm_density_divided_by_mean, geom, nStep, "rhodm", -5.0, 5.0, 200);
        total_time = ParallelDescriptor::second() - time1;
        if (ParallelDescriptor::IOProcessor())
            std::cout << std::setw(15) << " done. (" << total_time << " sec)" << std::endl;

        if (ParallelDescriptor::IOProcessor())
          std::cout << std::setfill('=') << std::setw(46) << " Calculating rho_dm P(k) ... " << std::flush;
        time1 = ParallelDescriptor::second();
        #warning TODO: make dm_density power spectrum knobs run-time parameters
        ps3d (dm_density, geom, nStep, "rhodm", true);
        total_time = ParallelDescriptor::second() - time1;
        if (ParallelDescriptor::IOProcessor())
            std::cout << std::setw(15) << " done. (" << total_time << " sec)" << std::endl;
    }

    {
        if (ParallelDescriptor::IOProcessor())
            std::cout << std::endl << std::setfill('=') << std::setw(46) << " Loading rho_m data ... " << std::flush;
        time1 = ParallelDescriptor::second();
        MultiFab rho_m (ba1, 1, 0);
        const int density_num_ghosts = density.nGrow();
        const int dm_density_num_ghosts = dm_density.nGrow();
        const int rho_m_num_ghosts = rho_m.nGrow();
        const Real mean_dm_density = dm_density.norm1() / Real(num_cells);
        for (MFIter mfi(dm_density); mfi.isValid(); ++mfi) {
          const Box& bx = mfi.validbox();
          BL_FORT_PROC_CALL(CALC_RHO_M, calc_rho_m) (
            (density)[mfi].dataPtr(),
            (dm_density)[mfi].dataPtr(),
            bx.loVect(),
            bx.hiVect(),
            &density_num_ghosts,
            &dm_density_num_ghosts,
            &rho_m_num_ghosts,
            &omega_b,
            &omega_m,
            &mean_density,
            &mean_dm_density,
            (rho_m)[mfi].dataPtr());
        }
        total_time = ParallelDescriptor::second() - time1;
        if (ParallelDescriptor::IOProcessor())
            std::cout << std::setw(15) << " done. (" << total_time << " sec)" << std::endl;

        if (ParallelDescriptor::IOProcessor())
          std::cout << std::setfill('=') << std::setw(46) << " Calculating rho_m PDF ... " << std::flush;
        time1 = ParallelDescriptor::second();
        #warning TODO: make rho_m PDF knobs run-time parameters
        pdf (rho_m, geom, nStep, "rhom", -3.0, 5.0, 200);
        total_time = ParallelDescriptor::second() - time1;
        if (ParallelDescriptor::IOProcessor())
            std::cout << std::setw(15) << " done. (" << total_time << " sec)" << std::endl;
    }

    {
        if (ParallelDescriptor::IOProcessor())
            std::cout << std::endl << std::setfill('=') << std::setw(46) << " Loading/deconvolving rho_dm data ... " << std::flush;
        time1 = ParallelDescriptor::second();

        // First query FFTW to get the slab size that it prefers for the domain
        // decomposition. We'll then construct boxes based on this distribution.
        BoxList bl_fft;
        intptr_t alloc_local;
        MakeFFTWBoxes(geom, bl_fft, alloc_local);
        BoxArray ba_fft(bl_fft);

        // Now that we have the domain decomposition that FFTW wants to use, we
        // set up a DM which places at most ONE Box on every process (which
        // FFTW implicitly assumes we will do; it has no notion of distributed
        // iterators). Presumably any of the supported distribution mapping
        // strategies in BoxLib will place one Box per process if there are
        // exactly as many Boxes as processors. But I don't know this for sure,
        // so we enforce it by hand.
        //
        // TODO: figure out if the DistributionMapping strategies try to put
        // more than 1 Box on any process if n_procs > n_boxes.
        const int nprocs = ParallelDescriptor::NProcs();
        Array<int> dm_one_box_per_process_vector(bl_fft.size()+1);
        for (unsigned int i = 0; i < bl_fft.size(); ++i)
            dm_one_box_per_process_vector[i] = i;
        dm_one_box_per_process_vector[bl_fft.size()] = ParallelDescriptor::MyProc();
        DistributionMapping dm_one_box_per_process(dm_one_box_per_process_vector);

        MultiFab rho_dm_fft(ba_fft, 1, 0, dm_one_box_per_process);
        rho_dm_fft.copy(dm_density);

        const Box& domain = geom.Domain();
        const IntVect& domain_size_IV = domain.size();
        const int* const domain_size_int = domain_size_IV.getVect();
        const Real dx = (geom.ProbHi(0) - geom.ProbLo(0)) / Real(domain.length(0));

        DistributionMapping dm(ba1, ParallelDescriptor::NProcs());

        // Get list of processes that own Boxes.
        int I_have_boxes = 0;
        for (MFIter mfi(rho_dm_fft); mfi.isValid(); ++mfi) {
          I_have_boxes++;
          break;
        }
        Array<int> proc_box_count(ParallelDescriptor::NProcs());
        MPI_Allgather(&I_have_boxes, 1, MPI_INT, &proc_box_count[0], 1, MPI_INT, ParallelDescriptor::Communicator());
        Array<int> procs_with_boxes;
        for (unsigned int i = 0; i < ParallelDescriptor::NProcs(); ++i) {
          if (proc_box_count[i] > 0) procs_with_boxes.push_back(i);
        }
        MPI_Group parent_group;
        MPI_Comm_group(ParallelDescriptor::Communicator(), &parent_group);
        MPI_Group fftw_group;
        MPI_Group_incl(parent_group, procs_with_boxes.size(), &procs_with_boxes[0], &fftw_group);
        MPI_Comm fftw_comm;
        MPI_Comm_create(ParallelDescriptor::Communicator(), fftw_group, &fftw_comm);

        const int comm = MPI_Comm_c2f(fftw_comm);

        // Separate MultiFabs to save real and imaginary components of complex
        // output from DFT. The Fortran implementation in BoxLib of MultiFabs
        // supports complex types but apparently the C++ implementation does not.
        MultiFab rho_dm_fft_out_real(ba_fft, 1, 0, dm_one_box_per_process);
        MultiFab rho_dm_fft_out_imag(ba_fft, 1, 0, dm_one_box_per_process);

        // Do the DFT. Because every process has at most one Box (or zero),
        // each process should iterate at most one time through this loop.
        for (MFIter mfi(rho_dm_fft); mfi.isValid(); ++mfi) {
            const Box bx = mfi.validbox();
            BL_FORT_PROC_CALL(FFT_3D, fft_3d) (
                (rho_dm_fft)[mfi].dataPtr(),
                bx.loVect(),
                bx.hiVect(),
                domain_size_int,
                &dx,
                &comm,
                &alloc_local,
                (rho_dm_fft_out_real)[mfi].dataPtr(),
                (rho_dm_fft_out_imag)[mfi].dataPtr());
        }

        // Do CIC convolution on k-space dark matter data
        const Box& problem_domain = geom.Domain();
        Array<int> domain_size(3);
        for (unsigned int i = 0; i < 3; ++i) domain_size[i] = problem_domain.length(i);
        for (MFIter mfi(rho_dm_fft_out_real); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();
            BL_FORT_PROC_CALL(CIC_DECONVOLVE, cic_deconvolve) (
                (rho_dm_fft_out_real)[mfi].dataPtr(),
                (rho_dm_fft_out_imag)[mfi].dataPtr(),
                bx.loVect(),
                bx.hiVect(),
                &domain_size[0]);
        }

        // Transform the CIC-deconvolved k-space data back to real space. The FTed
        // MultiFabs already have the correct Box decomposition that FFTW needs, so
        // we don't need to regrid.

        MultiFab rho_dm_cic_deconvolved_real(ba_fft, 1, 0, dm_one_box_per_process);
        MultiFab rho_dm_cic_deconvolved_imag(ba_fft, 1, 0, dm_one_box_per_process);

        for (MFIter mfi(rho_dm_fft_out_real); mfi.isValid(); ++mfi) {
            const Box bx = mfi.validbox();
            BL_FORT_PROC_CALL(FFT_3D_BACKWARD, fft_3d_backward) (
                (rho_dm_fft_out_real)[mfi].dataPtr(),
                (rho_dm_fft_out_imag)[mfi].dataPtr(),
                bx.loVect(),
                bx.hiVect(),
                domain_size_int,
                &dx,
                &comm,
                &alloc_local,
                (rho_dm_cic_deconvolved_real)[mfi].dataPtr(),
                (rho_dm_cic_deconvolved_imag)[mfi].dataPtr());
        }

        MPI_Group_free(&fftw_group);
        MPI_Group_free(&parent_group);
        if (MPI_COMM_NULL != fftw_comm) MPI_Comm_free(&fftw_comm);

        // We don't need the dumb slab decomposition anymore now that the DFT is
        // done. So re-grid onto regular boxes.

        MultiFab rho_dm_cic_deconvolved_real_regular_ba(ba1, 1, 0, dm);
        MultiFab rho_dm_cic_deconvolved_imag_regular_ba(ba1, 1, 0, dm);

        rho_dm_cic_deconvolved_real_regular_ba.copy(rho_dm_cic_deconvolved_real);
        rho_dm_cic_deconvolved_imag_regular_ba.copy(rho_dm_cic_deconvolved_imag);

        // Set matter overdensity using cosmological parameters. First get rho_m in
        // terms of mean units.

        const int density_num_ghosts = density.nGrow();
        const int rho_dm_num_ghosts = rho_dm_cic_deconvolved_real_regular_ba.nGrow();
        const long num_cells = problem_domain.numPts();
        const Real mean_density = density.norm1() / Real(num_cells);
        const Real mean_rho_dm = dm_density.norm1() / Real(num_cells);
        MultiFab rho_m(ba1, 1, 0);
        const int rho_m_num_ghosts = rho_m.nGrow();
        for (MFIter mfi(density); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();
            BL_FORT_PROC_CALL(CALC_RHO_M, calc_rho_m) (
              (density)[mfi].dataPtr(),
              (rho_dm_cic_deconvolved_real_regular_ba)[mfi].dataPtr(),
              bx.loVect(),
              bx.hiVect(),
              &density_num_ghosts,
              &rho_dm_num_ghosts,
              &rho_m_num_ghosts,
              &omega_b,
              &omega_m,
              &mean_density,
              &mean_rho_dm,
              (rho_m)[mfi].dataPtr());
        }
        // Now subtract 1 to get overdensity.
        rho_m.plus(-1.0, 0, 1);
        total_time = ParallelDescriptor::second() - time1;
        if (ParallelDescriptor::IOProcessor())
            std::cout << std::setw(15) << " done. (" << total_time << " sec)" << std::endl;

        if (ParallelDescriptor::IOProcessor())
          std::cout << std::setfill('=') << std::setw(46) << " Calculating rho_m P(k)... " << std::flush;
        time1 = ParallelDescriptor::second();
        #warning TODO: make rho_m power spectrum knobs run-time parameters
        ps3d (rho_m, geom, nStep, "rhom", false, true);
        total_time = ParallelDescriptor::second() - time1;
        if (ParallelDescriptor::IOProcessor())
            std::cout << std::setw(15) << " done. (" << total_time << " sec)" << std::endl;
    }

    {
        if (ParallelDescriptor::IOProcessor())
            std::cout << std::endl << std::setfill('=') << std::setw(46) << " Loading abs_v data ... " << std::flush;
        time1 = ParallelDescriptor::second();

        MultiFab abs_v(xmom.boxArray(), 1, 0);
        const int mom_num_ghosts = xmom.nGrow();
        const int density_num_ghosts = density.nGrow();

        for (MFIter mfi(xmom); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();
            BL_FORT_PROC_CALL(CALC_ABS_V, calc_abs_v) (
              (xmom)[mfi].dataPtr(),
              (ymom)[mfi].dataPtr(),
              (zmom)[mfi].dataPtr(),
              &mom_num_ghosts,
              (density)[mfi].dataPtr(),
              &density_num_ghosts,
              bx.loVect(),
              bx.hiVect(),
              (abs_v)[mfi].dataPtr());
        }
        total_time = ParallelDescriptor::second() - time1;
        if (ParallelDescriptor::IOProcessor())
            std::cout << std::setw(15) << " done. (" << total_time << " sec)" << std::endl;

        if (ParallelDescriptor::IOProcessor())
          std::cout << std::setfill('=') << std::setw(46) << " Calculating abs_v PDF ... " << std::flush;
        time1 = ParallelDescriptor::second();
        #warning TODO: make abs_v PDF knobs run-time parameters
        pdf (abs_v, geom, nStep, "velmag", -2.0, 4.0, 200);
        total_time = ParallelDescriptor::second() - time1;
        if (ParallelDescriptor::IOProcessor())
            std::cout << std::setw(15) << " done. (" << total_time << " sec)" << std::endl;

        if (ParallelDescriptor::IOProcessor())
          std::cout << std::setfill('=') << std::setw(46) << " Calculating abs_v P(k)... " << std::flush;
        time1 = ParallelDescriptor::second();
        #warning TODO: make abs_v power spectrum knobs run-time parameters
        ps3d(abs_v, geom, nStep, "velmag");
        total_time = ParallelDescriptor::second() - time1;
        if (ParallelDescriptor::IOProcessor())
            std::cout << std::setw(15) << " done. (" << total_time << " sec)" << std::endl;
    }

    {
        if (ParallelDescriptor::IOProcessor())
            std::cout << std::endl << std::setfill('=') << std::setw(46) << " Loading abs_vz data ... " << std::flush;
        time1 = ParallelDescriptor::second();
        MultiFab abs_vz(zmom.boxArray(), 1, 0);
        const int zmom_num_ghosts = zmom.nGrow();
        const int density_num_ghosts = density.nGrow();
        const int vz_num_ghosts = abs_vz.nGrow();
        for (MFIter mfi(zmom); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();
            BL_FORT_PROC_CALL(CALC_ABS_VZ, calc_abs_vz) (
                (zmom)[mfi].dataPtr(),
                (density)[mfi].dataPtr(),
                &zmom_num_ghosts,
                &density_num_ghosts,
                &vz_num_ghosts,
                bx.loVect(),
                bx.hiVect(),
                (abs_vz)[mfi].dataPtr());
        }
        total_time = ParallelDescriptor::second() - time1;
        if (ParallelDescriptor::IOProcessor())
            std::cout << std::setw(15) << " done. (" << total_time << " sec)" << std::endl;

        if (ParallelDescriptor::IOProcessor())
          std::cout << std::setfill('=') << std::setw(46) << " Calculating abs_vz PDF ... " << std::flush;
        time1 = ParallelDescriptor::second();
        #warning TODO: make abs_vz PDF knobs run-time parameters
        pdf (abs_vz, geom, nStep, "vz", -2.0, 4.0, 200);
        total_time = ParallelDescriptor::second() - time1;
        if (ParallelDescriptor::IOProcessor())
            std::cout << std::setw(15) << " done. (" << total_time << " sec)" << std::endl;

        if (ParallelDescriptor::IOProcessor())
          std::cout << std::setfill('=') << std::setw(46) << " Calculating abs_vz P(k) ... " << std::flush;
        time1 = ParallelDescriptor::second();
        #warning TODO: make abs_vz power spectrum knobs run-time parameters
        ps3d(abs_vz, geom, nStep, "vz");
        total_time = ParallelDescriptor::second() - time1;
        if (ParallelDescriptor::IOProcessor())
            std::cout << std::setw(15) << " done. (" << total_time << " sec)" << std::endl;
    }

    {
        if (ParallelDescriptor::IOProcessor())
            std::cout << std::endl << std::setfill('=') << std::setw(46) << " Loading n_hi data ... " << std::flush;
        time1 = ParallelDescriptor::second();
        MultiFab n_hi(density.boxArray(), 1, 0);
        const int state_num_ghosts = density.nGrow();
        const int n_hi_num_ghosts = n_hi.nGrow();
        for (MFIter mfi(n_hi); mfi.isValid(); ++mfi) {
          const Box& bx = mfi.validbox();
          BL_FORT_PROC_CALL(CALC_N_HI, calc_n_hi) (
            &z,
            (density)[mfi].dataPtr(),
            (e_int)[mfi].dataPtr(),
            &mean_density,
            &omega_b,
            &H_0,
            bx.loVect(),
            bx.hiVect(),
            &state_num_ghosts,
            &n_hi_num_ghosts,
            (n_hi)[mfi].dataPtr());
        }
        total_time = ParallelDescriptor::second() - time1;
        if (ParallelDescriptor::IOProcessor())
            std::cout << std::setw(15) << " done. (" << total_time << " sec)" << std::endl;

        if (ParallelDescriptor::IOProcessor())
          std::cout << std::setfill('=') << std::setw(46) << " Calculating n_hi PDF ... " << std::flush;
        time1 = ParallelDescriptor::second();
        #warning TODO: make n_hi PDF knobs run-time parameters
        pdf (n_hi, geom, nStep, "nhi", -13.0, -5.0, 200);
        total_time = ParallelDescriptor::second() - time1;
        if (ParallelDescriptor::IOProcessor())
            std::cout << std::setw(15) << " done. (" << total_time << " sec)" << std::endl;

        if (ParallelDescriptor::IOProcessor())
          std::cout << std::setfill('=') << std::setw(46) << " Calculating n_hi P(k) ... " << std::flush;
        time1 = ParallelDescriptor::second();
        #warning TODO: make n_hi power spectrum knobs run-time parameters
        ps3d (n_hi, geom, nStep, "nhi");
        total_time = ParallelDescriptor::second() - time1;
        if (ParallelDescriptor::IOProcessor())
            std::cout << std::setw(15) << " done. (" << total_time << " sec)" << std::endl;
    }
}
