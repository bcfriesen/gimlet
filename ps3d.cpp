#include <cmath>
#include <iomanip>
#include <fstream>
#include <sstream>

#include <DistributionMapping.H>
#include <MultiFab.H>
#include <Geometry.H>
#include <BLFort.H>

#include <MakeFFTWBoxes.H>

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

BL_FORT_PROC_DECL(CALC_PS3D, calc_ps3d) (
    const Real* mf_fft_out_real,
    const Real* mf_fft_out_imag,
    const int* lo,
    const int* hi,
    const int* num_ghosts,
    const int* num_bins,
    const Real* k_bin_edges,
    const Real* domain_length,
    const int* domain_grid_length,
    const int* k_bin_count,
    const Real* k_bin_power_weighted_k_sum,
    const Real* k_bin_power_sum);

BL_FORT_PROC_DECL(CIC_DECONVOLVE, cic_deconvolve) (
    const Real* dm_density_real,
    const Real* dm_density_imag,
    const int* lo,
    const int* hi,
    const int* domain_size);

void ps3d (const MultiFab& mf,
           const Geometry &geom,
           const unsigned int nStep,
           const std::string field_name,
           bool CIC_deconvolve,
           bool already_in_overdensity_units)
{
    /*
     * We can do pencil FFTs on single processors (i.e., without MPI) because
     * the pencil memory footprint is small. However this is not true for the
     * grid-wide 3-D FFTs. For these we must use the distributed memory MPI
     * version of FFTW.
     */

    // First query FFTW to get the slab size that it prefers for the domain
    // decomposition. We'll then construct boxes based on this distribution.
    BoxList bl_fft;
    intptr_t alloc_local;
    MakeFFTWBoxes(geom, bl_fft, alloc_local);

    // Now that we have the domain decomposition that FFTW wants to use, we set
    // up a DM which places exactly one Box on every process (which FFTW
    // implicitly assumes we will do; it has no notion of distributed
    // iterators). Presumably any of the supported distribution mapping
    // strategies in BoxLib will place one Box per process if there are exactly
    // as many Boxes as processors. But I don't know this for sure, so we
    // enforce it by hand.
    Array<int> dm_one_box_per_process_vector(bl_fft.size()+1);
    for (unsigned int i = 0; i < bl_fft.size(); ++i)
        dm_one_box_per_process_vector[i] = i;
    dm_one_box_per_process_vector[bl_fft.size()] = ParallelDescriptor::MyProc();
    DistributionMapping dm_one_box_per_process(dm_one_box_per_process_vector);

    BoxArray ba_fft(bl_fft);
    MultiFab mf_fft(ba_fft, 1, 0, dm_one_box_per_process);

    mf_fft.copy(mf);

    const Box& domain = geom.Domain();
    const long num_cells = domain.numPts();

    MultiFab overmf(mf_fft.boxArray(), 1, 0, dm_one_box_per_process);
    overmf.copy(mf_fft);

    if (!already_in_overdensity_units) {
        const Real mean_value = mf.norm1() / Real(num_cells);
        overmf.mult(1.0/mean_value);
        overmf.plus(-1.0, 0, 1);
    }

    const IntVect& domain_size_IV = domain.size();
    const int* const domain_size_int = domain_size_IV.getVect();
    const Real dx = (geom.ProbHi(0) - geom.ProbLo(0)) / Real(domain.length(0));

    // The following section is devoted to the case that there are more
    // processes than there are slabs for FFTW to distribute (i.e., n_procs >
    // grid_nz). Since the calls to MPI FFTW are collective, but not all
    // processes will make the call (since some won't have any slabs to work
    // on), we confine the "active" processes to a subcommunicator during the
    // DFT.

    // Get list of processes that have a slab to work on.
    int num_my_boxes = 0;
    for (MFIter mfi(overmf); mfi.isValid(); ++mfi) {
      num_my_boxes++;
      break;
    }
    Array<int> proc_box_count(ParallelDescriptor::NProcs());
    MPI_Allgather(&num_my_boxes, 1, MPI_INT, &proc_box_count[0], 1, MPI_INT, ParallelDescriptor::Communicator());
    Array<int> procs_with_boxes;
    for (unsigned int i = 0; i < ParallelDescriptor::NProcs(); ++i) {
      if (proc_box_count[i] > 0) procs_with_boxes.push_back(i);
    }
    // Build subgroup and subcommunicator for those processes which will
    // participate in the DFT.
    MPI_Group parent_group;
    MPI_Comm_group(ParallelDescriptor::Communicator(), &parent_group);
    MPI_Group fftw_group;
    MPI_Group_incl(parent_group, procs_with_boxes.size(), &procs_with_boxes[0], &fftw_group);
    MPI_Comm fftw_comm;
    MPI_Comm_create(ParallelDescriptor::Communicator(), fftw_group, &fftw_comm);

    const int comm = MPI_Comm_c2f(fftw_comm);

    // Separate MultiFabs to save real and imaginary components of complex
    // output from DFT.
    MultiFab mf_fft_out_real(ba_fft, 1, 0, dm_one_box_per_process);
    MultiFab mf_fft_out_imag(ba_fft, 1, 0, dm_one_box_per_process);

    // Do the DFT. Every process has either 1 or 0 boxes to work on, so even
    // though we use a generic MFIter, it will only iterate one time (or zero).
    for (MFIter mfi(overmf); mfi.isValid(); ++mfi) {
        const Box bx = mfi.validbox();
        BL_FORT_PROC_CALL(FFT_3D, fft_3d) (
            (overmf)[mfi].dataPtr(),
            bx.loVect(),
            bx.hiVect(),
            domain_size_int,
            &dx,
            &comm,
            &alloc_local,
            (mf_fft_out_real)[mfi].dataPtr(),
            (mf_fft_out_imag)[mfi].dataPtr());
    }

    MPI_Group_free(&fftw_group);
    MPI_Group_free(&parent_group);
    if (MPI_COMM_NULL != fftw_comm) MPI_Comm_free(&fftw_comm);

    if (CIC_deconvolve) {
        // Do CIC deconvolution on k-space dark matter data.
        const Box& problem_domain = geom.Domain();
        Array<int> domain_size(3);
        for (unsigned int i = 0; i < 3; ++i) domain_size[i] = problem_domain.length(i);
        for (MFIter mfi(mf_fft_out_real); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();
            BL_FORT_PROC_CALL(CIC_DECONVOLVE, cic_deconvolve) (
                (mf_fft_out_real)[mfi].dataPtr(),
                (mf_fft_out_imag)[mfi].dataPtr(),
                bx.loVect(),
                bx.hiVect(),
                &domain_size[0]);
        }
    }

    // Now calculate 3-D power spectrum P(k).

    const int n_min = std::min(domain.length(0), std::min(domain.length(1), domain.length(2)));
    const int i_nyq = n_min / 2;
    const Real pi = 3.1415927;
    const Real k_fund = (2.0 * pi) / geom.ProbLength(0);
    const Real k_nyquist = k_fund * Real(n_min)/2.0;

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
        } else {
            ki = i - num_lin_bins;
            k_bin_edges[i] = pow(10.0, log_k_lin_edge + dlogk * ki);
        }
    }

    if (k_min < 0.0 || k_max > k_nyquist)
      BoxLib::Error("Bad k bin edge values.");

    const int num_ghosts = mf_fft_out_real.nGrow();
    const Real domain_length = geom.ProbLength(0);

    std::vector<int>  k_bin_count                (k_num_bins);
    std::vector<Real> k_bin_power_weighted_k_sum (k_num_bins);
    std::vector<Real> k_bin_power_sum            (k_num_bins);

    for (unsigned int i = 0; i < k_num_bins; ++i) {
        k_bin_count[i] = 0;
        k_bin_power_weighted_k_sum[i] = 0.0;
        k_bin_power_sum[i] = 0.0;
    }

    const int domain_grid_length = domain.length(0);

    for (MFIter mfi(mf_fft_out_real); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();

        std::vector<int>  k_bin_count_per_box                (k_num_bins);
        std::vector<Real> k_bin_power_weighted_k_sum_per_box (k_num_bins);
        std::vector<Real> k_bin_power_sum_per_box            (k_num_bins);

        BL_FORT_PROC_CALL(CALC_PS3D, calc_ps3d) (
            (mf_fft_out_real)[mfi].dataPtr(),
            (mf_fft_out_imag)[mfi].dataPtr(),
            bx.loVect(),
            bx.hiVect(),
            &num_ghosts,
            &k_num_bins,
            &k_bin_edges[0],
            &domain_length,
            &domain_grid_length,
            &k_bin_count[0],
            &k_bin_power_weighted_k_sum[0],
            &k_bin_power_sum[0]);

        for (unsigned int i = 0; i < k_num_bins; ++i) {
             k_bin_count[i] += k_bin_count_per_box[i];
             k_bin_power_weighted_k_sum[i] += k_bin_power_weighted_k_sum_per_box[i];
             k_bin_power_sum[i] += k_bin_power_sum_per_box[i];
        }
    }

    ParallelDescriptor::ReduceIntSum (&k_bin_count[0], k_num_bins);
    ParallelDescriptor::ReduceRealSum(&k_bin_power_weighted_k_sum[0], k_num_bins);
    ParallelDescriptor::ReduceRealSum(&k_bin_power_sum[0], k_num_bins);

    if (ParallelDescriptor::IOProcessor()) {
        std::ofstream power_spectrum_file;
        std::stringstream filename;
        filename << field_name << "_power_spectrum_pdf_nstep_" << std::setfill('0') << std::setw(5) << nStep;
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
    return;
}
