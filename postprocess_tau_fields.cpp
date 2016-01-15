#include <fstream>
#include <sstream>
#include <iomanip>

#include <MultiFab.H>
#include <Geometry.H>
#include <BLFort.H>

#include <MakePencilBoxes.H>

BL_FORT_PROC_DECL(CALC_PDF, calc_pdf) (
        Real       *tau,
        const int  *lo,
        const int  *hi,
        const int  *num_ghosts,
        const int  *pdf_num_bins,
        const Real *bin_edges,
        int        *bin_count,
        Real       *bin_x_sum);

BL_FORT_PROC_DECL(CALC_FLUX, calc_flux) (
        Real       *tau,
        const int  *lo,
        const int  *hi,
        const int  *num_ghosts,
        const Real *A,
        Real       *flux);

BL_FORT_PROC_DECL(CALC_OVERDENSITY, calc_overdensity) (
        Real       *density,
        const int  *lo,
        const int  *hi,
        const int  *num_ghosts,
        const Real *mean_density_inv,
        Real       *overdensity_real);

BL_FORT_PROC_DECL(CALC_PENCIL_FFT, calc_pencil_fft) (
        Real       *overdensity,
        const int  *lo,
        const int  *hi,
        const int  *num_ghosts,
        const int  *dir,
        const Real *domain_length,
        const Real *dx,
        const int  *num_bins,
        const Real *k_bin_edges,
        Real       *overdensity_fft_real,
        Real       *overdensity_fft_imag);

BL_FORT_PROC_DECL(CALC_PS1D, calc_ps1d) (
        Real       *overdensity_fft_real,
        Real       *overdensity_fft_imag,
        const int  *lo,
        const int  *hi,
        const int  *ng,
        const int  *dir,
        const int  *num_bins,
        const Real *k_bin_edges,
        const Real *domain_length,
        int        *k_bin_count,
        Real       *k_bin_power_weighted_k_sum,
        Real       *k_bin_power_sum);

BL_FORT_PROC_DECL(CALC_K_MU_PS1D, calc_k_mu_ps1d) (
        Real       *overdensity_fft_real,
        Real       *overdensity_fft_imag,
        const int  *lo,
        const int  *hi,
        const int  *ng,
        const int  *dir,
        const Real *domain_length,
        const int  *k_num_bins,
        const Real *k_bin_edges,
        const int  *mu_num_bins,
        const Real *mu_bin_edges,
        int        *k_bin_count,
        Real       *k_bin_power_weighted_k_sum,
        Real       *k_bin_power_weighted_mu_sum,
        Real       *k_bin_power_sum);

void postprocess_tau_fields(MultiFab& tau, const Geometry &geom, const int dir, const int nStep)
{

    static bool first_time_writing_mean_flux_file_this_timestep = true;

#warning TODO: make PDF knobs run-time parameters
    const Real lt_min = -2.0;
    const Real lt_max =  3.0;
    const int tpdf_num_bins = 200;
    const int fpdf_num_bins = 50;

    const int num_ghosts = 0;

    // These are counts per *process*, not per *box*. We'll AllReduce() them
    // after iterating over all boxes.
    std::vector<int>  bin_count;
    std::vector<Real> bin_x_sum;
    std::vector<Real> bin_edges;

    bin_count.resize(tpdf_num_bins);
    bin_x_sum.resize(tpdf_num_bins);
    for (unsigned int i = 0; i < tpdf_num_bins; ++i) {
        bin_count[i] = 0;
        bin_x_sum[i] = 0.0;
    }

    bin_edges.resize(tpdf_num_bins+1);
    Real bin_width = (lt_max - lt_min) / Real(tpdf_num_bins);
    for (unsigned int i = 0; i <= tpdf_num_bins; ++i) {
        bin_edges[i] = std::pow(10.0, lt_min + Real(i)*bin_width);
    }

    if (ParallelDescriptor::IOProcessor())
      std::cout << std::setfill('=') << std::setw(46) << " Calculating tau PDF ... " << std::flush;

    Real time1 = ParallelDescriptor::second();
    for ( MFIter mfi(tau); mfi.isValid(); ++mfi ) {
        const Box& bx = mfi.validbox();

        // These are counts per *box*. We'll glob them all together after
        // iterating over all boxes.
        std::vector<int>  bin_count_per_box(tpdf_num_bins);
        std::vector<Real> bin_x_sum_per_box(tpdf_num_bins);

        BL_FORT_PROC_CALL(CALC_PDF, calc_pdf) (
                (tau)[mfi].dataPtr(),
                bx.loVect(),
                bx.hiVect(),
                &num_ghosts,
                &tpdf_num_bins,
                &bin_edges[0],
                &bin_count_per_box[0],
                &bin_x_sum_per_box[0]);

        // Accumulate per-box results into per-process arrays
        for (unsigned int i = 0; i < tpdf_num_bins; ++i) {
             bin_count[i] += bin_count_per_box[i];
             bin_x_sum[i] += bin_x_sum_per_box[i];
        }
    }

    ParallelDescriptor::ReduceIntSum (&bin_count[0], tpdf_num_bins);
    ParallelDescriptor::ReduceRealSum(&bin_x_sum[0], tpdf_num_bins);

    // Normalization
    Real sum_ni_dxi = 0.0;
    for (unsigned int i = 0; i < tpdf_num_bins; ++i) {
        Real dx_i = bin_edges[i+1] - bin_edges[i];
        sum_ni_dxi += bin_count[i] * dx_i;
    }

    // Normalization constant for N_i -> P_i
    Real pnorm = 1.0 / sum_ni_dxi;

    std::vector<Real> bin_x(tpdf_num_bins);
    std::vector<Real> bin_p(tpdf_num_bins);

    for (unsigned int i = 0; i < tpdf_num_bins; ++i) {
        if (bin_count[i] > 0) {
            bin_x[i] = bin_x_sum[i] / Real(bin_count[i]);
            bin_p[i] = pnorm * Real(bin_count[i]);
        }
    }
    Real total_time = ParallelDescriptor::second() - time1;
    if (ParallelDescriptor::IOProcessor())
      std::cout << std::setw(15) << " done. (" << total_time << " sec)" << std::endl;

    if (ParallelDescriptor::IOProcessor()) {

        std::ofstream tau_pdf_output_file;
        std::stringstream filename;
        filename << "tau_pdf_nstep_" << std::setfill('0') << std::setw(5) << nStep << "_dir_" << dir;
        tau_pdf_output_file.open(filename.str().c_str());
        tau_pdf_output_file << std::scientific;
        for (unsigned int i = 0; i < tpdf_num_bins; ++i) {
            tau_pdf_output_file << std::setw(25) << bin_edges[i] << std::setw(15) << bin_count[i] << std::setw(15) << bin_x[i] << std::setw(15) << bin_p[i] << std::endl;
        }
        tau_pdf_output_file.close();
    }

    // Compute flux fields. F ~ exp(-tau)
    const BoxArray& ba = tau.boxArray();
    MultiFab flux(ba, 1, num_ghosts);

#warning TODO: MAKE FLUX SCALE FACTOR A KNOB
#warning TODO: ADD OTHER TWO MODES FOR FLUX CALCULATION

    const Real A = 1.0;

    if (ParallelDescriptor::IOProcessor())
      std::cout << std::setfill('=') << std::setw(46) << " Calculating fluxes ... " << std::flush;

    time1 = ParallelDescriptor::second();
    for ( MFIter mfi(tau); mfi.isValid(); ++mfi ) {
        const Box& bx = mfi.validbox();
        BL_FORT_PROC_CALL(CALC_FLUX, calc_flux) (
                (tau)[mfi].dataPtr(),
                bx.loVect(),
                bx.hiVect(),
                &num_ghosts,
                &A,
                (flux)[mfi].dataPtr());
    }
    total_time = ParallelDescriptor::second() - time1;
    if (ParallelDescriptor::IOProcessor())
      std::cout << std::setw(15) << " done. (" << total_time << " sec)" << std::endl;

    Real mean_flux;

    if (ParallelDescriptor::IOProcessor())
      std::cout << std::setfill('=') << std::setw(46) << " Calculating mean flux ... " << std::flush;

    time1 = ParallelDescriptor::second();
    mean_flux = flux.norm1() / Real(ba.numPts());

    total_time = ParallelDescriptor::second() - time1;
    if (ParallelDescriptor::IOProcessor())
      std::cout << std::setw(15) << " done. (" << total_time << " sec)" << std::endl;

    if (ParallelDescriptor::IOProcessor()) {
        std::ofstream mean_flux_output_file;
        std::stringstream filename;
        filename << "mean_fluxes_" << std::setfill('0') << std::setw(5) << nStep;
        if (first_time_writing_mean_flux_file_this_timestep) {
            mean_flux_output_file.open(filename.str().c_str());
            mean_flux_output_file << "#" << std::setw(14) << "direction" << std::setw(15) << "<F>" << std::endl;
            first_time_writing_mean_flux_file_this_timestep = false;
        } else {
            mean_flux_output_file.open(filename.str().c_str(), std::ios::app);
        }
        mean_flux_output_file << std::scientific;
        mean_flux_output_file << std::setw(15) << dir << std::setw(15) << mean_flux << std::endl;
        mean_flux_output_file.close();
    }
    // After we've written the mean fluxes in the z-direction, we'll move on to
    // the next timestep.
    if (dir == 2) first_time_writing_mean_flux_file_this_timestep = true;

    const Real flux_min = 0.0;
    const Real flux_max = 1.0;
    bin_width = (flux_max - flux_min) / Real(fpdf_num_bins);
    bin_edges.resize(fpdf_num_bins+1);
    for (unsigned int i = 0; i <= fpdf_num_bins; ++i) {
        bin_edges[i] = flux_min + Real(i)*bin_width;
    }

    bin_count.resize(fpdf_num_bins);
    bin_x_sum.resize(fpdf_num_bins);
    for (unsigned int i = 0; i < fpdf_num_bins; ++i) {
        bin_count[i] = 0;
        bin_x_sum[i] = 0.0;
    }

    if (ParallelDescriptor::IOProcessor())
      std::cout << std::setfill('=') << std::setw(46) << " Calculating flux PDF ... " << std::flush;

    time1 = ParallelDescriptor::second();
    for ( MFIter mfi(flux); mfi.isValid(); ++mfi ) {
        const Box& bx = mfi.validbox();

        // These are counts per *box*. We'll glob them all together after
        // iterating over all boxes.
        std::vector<int>  bin_count_per_box(fpdf_num_bins);
        std::vector<Real> bin_x_sum_per_box(fpdf_num_bins);

        BL_FORT_PROC_CALL(CALC_PDF, calc_pdf) (
                (flux)[mfi].dataPtr(),
                bx.loVect(),
                bx.hiVect(),
                &num_ghosts,
                &fpdf_num_bins,
                &bin_edges[0],
                &bin_count_per_box[0],
                &bin_x_sum_per_box[0]);

        for (unsigned int i = 0; i < fpdf_num_bins; ++i) {
             bin_count[i] += bin_count_per_box[i];
             bin_x_sum[i] += bin_x_sum_per_box[i];
        }
    }

    ParallelDescriptor::ReduceIntSum (&bin_count[0], fpdf_num_bins);
    ParallelDescriptor::ReduceRealSum(&bin_x_sum[0], fpdf_num_bins);

    // Normalization
    sum_ni_dxi = 0.0;
    for (unsigned int i = 0; i < fpdf_num_bins; ++i) {
        Real dx_i = bin_edges[i+1] - bin_edges[i];
        sum_ni_dxi += bin_count[i] * dx_i;
    }

    // Normalization constant for N_i -> P_i
    pnorm = 1.0 / sum_ni_dxi;

    for (unsigned int i = 0; i < fpdf_num_bins; ++i) {
        if (bin_count[i] > 0) {
            bin_x[i] = bin_x_sum[i] / Real(bin_count[i]);
            bin_p[i] = pnorm * Real(bin_count[i]);
        }
    }
    total_time = ParallelDescriptor::second() - time1;
    if (ParallelDescriptor::IOProcessor())
      std::cout << std::setw(15) << " done. (" << total_time << " sec)" << std::endl;

    if (ParallelDescriptor::IOProcessor()) {

        std::ofstream flux_pdf_output_file;
        std::stringstream filename;
        filename << "flux_pdf_nstep_" << std::setfill('0') << std::setw(5) << nStep << "_dir_" << dir;
        flux_pdf_output_file.open(filename.str().c_str());
        flux_pdf_output_file << std::scientific;
        for (unsigned int i = 0; i < fpdf_num_bins; ++i) {
            flux_pdf_output_file << std::setw(25) << bin_edges[i] << std::setw(15) << bin_count[i] << std::setw(15) << bin_x[i] << std::setw(15) << bin_p[i] << std::endl;
        }
        flux_pdf_output_file.close();
    }


    // Power spectra

    const Box& problem_domain = geom.Domain();
    const int n_min = std::min(problem_domain.length(0), std::min(problem_domain.length(1), problem_domain.length(2)));
    const int i_nyq = n_min / 2;
    const Real pi = 3.1415927;
    const Real k_fund = (2.0 * pi) / geom.ProbLength(0);

    // bin bounds
    const Real k_min = k_fund;
    const Real k_max = k_fund * Real(i_nyq);

    // Chose number of linear bins at low end.
    const int num_lin_bins = 10;
    const Real k_lin_edge = k_min + Real(num_lin_bins) * k_fund;

    // Figure out log spacing.
    // Choose dlogk such that the spacing between i = 10 and 11 looks like the
    // linear spacing would on log scale.
    const Real dlogk = log10( k_lin_edge + k_fund ) - log10( k_lin_edge );
    const int num_log_bins = ( log10(k_max) - log10(k_lin_edge) ) / dlogk;
    const Real log_k_lin_edge = log10(k_lin_edge);

    const int k_num_bins = num_lin_bins + num_log_bins;
    const int k_num_edges = k_num_bins + 1;

    std::vector<Real> k_bin_edges(k_num_edges);
    int i, ki;
    for (i = 0; i < k_num_edges; ++i) {
        if (i < num_lin_bins) {
            k_bin_edges[i] = k_min + k_fund * i;
        }
        else {
            ki = i - num_lin_bins;
            k_bin_edges[i] = pow(10.0, log_k_lin_edge + dlogk * ki);
        }
    }

    // Power spectra of 1-D rays

    const Real mean_flux_inv = 1.0 / mean_flux;

    MultiFab overflux(ba, 1, num_ghosts);

    for ( MFIter mfi(flux); mfi.isValid(); ++mfi ) {
        const Box& bx = mfi.validbox();

        BL_FORT_PROC_CALL(CALC_OVERDENSITY, calc_overdensity) (
                (flux)[mfi].dataPtr(),
                bx.loVect(),
                bx.hiVect(),
                &num_ghosts,
                &mean_flux_inv,
                (overflux)[mfi].dataPtr());
    }

    // Set up pencil boxes

    BoxList bl = MakePencilBoxes(geom, dir);
    BoxArray ba_pencil(bl);

    MultiFab overflux_pencils          (ba_pencil, 1, num_ghosts);
    MultiFab overflux_fft_pencils_real (ba_pencil, 1, num_ghosts);
    MultiFab overflux_fft_pencils_imag (ba_pencil, 1, num_ghosts);

    overflux_pencils.copy(overflux);

    const Real domain_length = geom.ProbLength(0);
    const Real dx = (geom.ProbHi(0) - geom.ProbLo(0)) / Real(problem_domain.length(0));

    if (ParallelDescriptor::IOProcessor())
      std::cout << std::setfill('=') << std::setw(46) << " Calculating pencil FFTs ... " << std::flush;

    time1 = ParallelDescriptor::second();
    for ( MFIter mfi(overflux_pencils); mfi.isValid(); ++mfi ) {
        const Box& bx = mfi.validbox();

        BL_FORT_PROC_CALL(CALC_PENCIL_FFT, calc_pencil_fft) (
                (overflux_pencils)[mfi].dataPtr(),
                bx.loVect(),
                bx.hiVect(),
                &num_ghosts,
                &dir,
                &domain_length,
                &dx,
                &k_num_bins,
                &k_bin_edges[0],
                (overflux_fft_pencils_real)[mfi].dataPtr(),
                (overflux_fft_pencils_imag)[mfi].dataPtr());
    }
    total_time = ParallelDescriptor::second() - time1;
    if (ParallelDescriptor::IOProcessor())
      std::cout << std::setw(15) << " done. (" << total_time << " sec)" << std::endl;

    std::vector<int>  k_bin_count                (k_num_bins);
    std::vector<Real> k_bin_power_weighted_k_sum (k_num_bins);
    std::vector<Real> k_bin_power_sum            (k_num_bins);

    for (unsigned int i = 0; i < k_num_bins; ++i) {
        k_bin_count[i] = 0;
        k_bin_power_weighted_k_sum[i] = 0.0;
        k_bin_power_sum[i] = 0.0;
    }

    if (ParallelDescriptor::IOProcessor())
      std::cout << std::setfill('=') << std::setw(46) << " Calculating 1D P(k) ... " << std::flush;

    time1 = ParallelDescriptor::second();
    for ( MFIter mfi(overflux_fft_pencils_real); mfi.isValid(); ++mfi ) {
        const Box& bx = mfi.validbox();

        std::vector<int>  k_bin_count_per_box                (k_num_bins);
        std::vector<Real> k_bin_power_weighted_k_sum_per_box (k_num_bins);
        std::vector<Real> k_bin_power_sum_per_box            (k_num_bins);

        BL_FORT_PROC_CALL(CALC_PS1D, calc_ps1d) (
                (overflux_fft_pencils_real)[mfi].dataPtr(),
                (overflux_fft_pencils_imag)[mfi].dataPtr(),
                bx.loVect(),
                bx.hiVect(),
                &num_ghosts,
                &dir,
                &k_num_bins,
                &k_bin_edges[0],
                &domain_length,
                &k_bin_count_per_box[0],
                &k_bin_power_weighted_k_sum_per_box[0],
                &k_bin_power_sum_per_box[0]);

        for (unsigned int i = 0; i < k_num_bins; ++i) {
             k_bin_count[i] += k_bin_count_per_box[i];
             k_bin_power_weighted_k_sum[i] += k_bin_power_weighted_k_sum_per_box[i];
             k_bin_power_sum[i] += k_bin_power_sum_per_box[i];
        }
    }

    ParallelDescriptor::ReduceIntSum (&k_bin_count[0], k_num_bins);
    ParallelDescriptor::ReduceRealSum(&k_bin_power_weighted_k_sum[0], k_num_bins);
    ParallelDescriptor::ReduceRealSum(&k_bin_power_sum[0], k_num_bins);

    total_time = ParallelDescriptor::second() - time1;
    if (ParallelDescriptor::IOProcessor())
      std::cout << std::setw(15) << " done. (" << total_time << " sec)" << std::endl;


    if (ParallelDescriptor::IOProcessor()) {
        std::ofstream power_spectrum_file;
        std::stringstream filename;
        filename << "flux_ps1d_nstep_" << std::setfill('0') << std::setw(5) << nStep << "_dir_" << dir;
        power_spectrum_file.open(filename.str().c_str());
        power_spectrum_file << std::scientific;
        for (unsigned int i = 0; i < k_num_bins; ++i) {
            power_spectrum_file << std::setw(25) << k_bin_edges[i]
                                << std::setw(15) << k_bin_count[i]
                                << std::setw(25) << k_bin_power_weighted_k_sum[i]/k_bin_power_sum[i]
                                << std::setw(15) << k_bin_power_sum[i]/Real(k_bin_count[i])
                                << std::endl;
        }
        power_spectrum_file.close();
    }


    const int mu_num_bins = 4;
    std::vector<Real> mu_bin_edges(mu_num_bins+1);
    const Real mu_max_edge = 1.0;
    const Real mu_min_edge = 0.0;
    bin_width = (mu_max_edge - mu_min_edge) / Real(mu_num_bins);
    for (int i = 0; i <= mu_num_bins; ++i) {
        mu_bin_edges[i] = mu_min_edge + Real(i)*bin_width;
    }
     const int num_k_and_mu_bins = mu_num_bins * k_num_bins;

    std::vector<Real> k_bin_power_weighted_mu_sum;
    k_bin_count.resize(num_k_and_mu_bins);
    k_bin_power_weighted_k_sum.resize(num_k_and_mu_bins);
    k_bin_power_weighted_mu_sum.resize(num_k_and_mu_bins);
    k_bin_power_sum.resize(num_k_and_mu_bins);

    for (unsigned int i = 0; i < num_k_and_mu_bins; ++i) {
        k_bin_count[i] = 0;
        k_bin_power_weighted_k_sum[i] = 0.0;
        k_bin_power_weighted_mu_sum[i] = 0.0;
        k_bin_power_sum[i] = 0.0;
    }

    if (ParallelDescriptor::IOProcessor())
      std::cout << std::setfill('=') << std::setw(46) << " Calculating P(k,mu) ... " << std::flush;

    time1 = ParallelDescriptor::second();
    for ( MFIter mfi(overflux_fft_pencils_real); mfi.isValid(); ++mfi ) {
        const Box& bx = mfi.validbox();

        std::vector<int>  k_bin_count_per_box                 (num_k_and_mu_bins);
        std::vector<Real> k_bin_power_weighted_k_sum_per_box  (num_k_and_mu_bins);
        std::vector<Real> k_bin_power_weighted_mu_sum_per_box (num_k_and_mu_bins);
        std::vector<Real> k_bin_power_sum_per_box             (num_k_and_mu_bins);

        BL_FORT_PROC_CALL(CALC_K_MU_PS1D, calc_k_mu_ps1d) (
                (overflux_fft_pencils_real)[mfi].dataPtr(),
                (overflux_fft_pencils_imag)[mfi].dataPtr(),
                bx.loVect(),
                bx.hiVect(),
                &num_ghosts,
                &dir,
                &domain_length,
                &k_num_bins,
                &k_bin_edges[0],
                &mu_num_bins,
                &mu_bin_edges[0],
                &k_bin_count_per_box[0],
                &k_bin_power_weighted_k_sum_per_box[0],
                &k_bin_power_weighted_mu_sum_per_box[0],
                &k_bin_power_sum_per_box[0]);

        for (unsigned int i = 0; i < num_k_and_mu_bins; ++i) {
             k_bin_count[i] += k_bin_count_per_box[i];
             k_bin_power_weighted_k_sum[i] += k_bin_power_weighted_k_sum_per_box[i];
             k_bin_power_weighted_mu_sum[i] += k_bin_power_weighted_mu_sum_per_box[i];
             k_bin_power_sum[i] += k_bin_power_sum_per_box[i];
        }
    }

    ParallelDescriptor::ReduceIntSum (&k_bin_count[0], num_k_and_mu_bins);
    ParallelDescriptor::ReduceRealSum(&k_bin_power_weighted_k_sum[0], num_k_and_mu_bins);
    ParallelDescriptor::ReduceRealSum(&k_bin_power_weighted_mu_sum[0], num_k_and_mu_bins);
    ParallelDescriptor::ReduceRealSum(&k_bin_power_sum[0], num_k_and_mu_bins);

    total_time = ParallelDescriptor::second() - time1;
    if (ParallelDescriptor::IOProcessor())
      std::cout << std::setw(15) << " done. (" << total_time << " sec)" << std::endl;

    if (ParallelDescriptor::IOProcessor()) {
        std::ofstream k_mu_power_spectrum_file;
        std::stringstream filename;
        filename << "flux_pkmu_nstep_" << std::setfill('0') << std::setw(5) << nStep << "_dir_" << dir;
        k_mu_power_spectrum_file.open(filename.str().c_str());
        k_mu_power_spectrum_file << std::scientific;
        for (unsigned int i = 0; i < k_num_bins; ++i) {
            for (unsigned int j = 0; j < mu_num_bins; ++j) {
                unsigned int k_mu_index = i*mu_num_bins + j;
                k_mu_power_spectrum_file << std::setw(25) << k_bin_edges[i]
                                    << std::setw(25) << mu_bin_edges[j]
                                    << std::setw(15) << k_bin_count[k_mu_index]
                                    << std::setw(25) << k_bin_power_weighted_k_sum[k_mu_index]/k_bin_power_sum[k_mu_index]
                                    << std::setw(25) << k_bin_power_weighted_mu_sum[k_mu_index]/k_bin_power_sum[k_mu_index]
                                    << std::setw(15) << k_bin_power_sum[k_mu_index]/Real(k_bin_count[k_mu_index])
                                    << std::endl;
            }
        }
        k_mu_power_spectrum_file.close();
    }

}
