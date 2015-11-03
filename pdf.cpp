#include <iomanip>
#include <fstream>
#include <sstream>

#include <MultiFab.H>
#include <Geometry.H>
#include <BLFort.H>

BL_FORT_PROC_DECL(CALC_PDF, calc_pdf) (
        const Real* mf,
        const int*  lo,
        const int*  hi,
        const int*  num_ghosts,
        const int*  pdf_num_bins,
        const Real* bin_edges,
        const int*  bin_count,
        const Real* bin_x_sum);

void pdf (const MultiFab& mf,
          const Geometry &geom,
          const unsigned int nStep,
          const std::string field_name,
          const Real log10_min,
          const Real log10_max,
          const int pdf_num_bins)
{

    const int mf_num_ghosts = mf.nGrow();

    // These are counts per *process*, not per *box*. We'll AllReduce() them
    // after iterating over all boxes.
    std::vector<int>  bin_count;
    std::vector<Real> bin_x_sum;
    std::vector<Real> bin_edges;

    bin_count.resize(pdf_num_bins);
    bin_x_sum.resize(pdf_num_bins);
    for (unsigned int i = 0; i < pdf_num_bins; ++i) {
        bin_count[i] = 0;
        bin_x_sum[i] = 0.0;
    }

    bin_edges.resize(pdf_num_bins+1);
    Real bin_width = (log10_max - log10_min) / Real(pdf_num_bins);
    for (unsigned int i = 0; i <= pdf_num_bins; ++i) {
        bin_edges[i] = std::pow(10.0, log10_min + Real(i)*bin_width);
    }

    for ( MFIter mfi(mf); mfi.isValid(); ++mfi ) {
        const Box& bx = mfi.validbox();

        // These are counts per *box*. We'll glob them all together after
        // iterating over all boxes.
        std::vector<int>  bin_count_per_box(pdf_num_bins);
        std::vector<Real> bin_x_sum_per_box(pdf_num_bins);

        BL_FORT_PROC_CALL(CALC_PDF, calc_pdf) (
                (mf)[mfi].dataPtr(),
                bx.loVect(),
                bx.hiVect(),
                &mf_num_ghosts,
                &pdf_num_bins,
                &bin_edges[0],
                &bin_count_per_box[0],
                &bin_x_sum_per_box[0]);

        // Accumulate per-box results into per-process arrays
        for (unsigned int i = 0; i < pdf_num_bins; ++i) {
             bin_count[i] += bin_count_per_box[i];
             bin_x_sum[i] += bin_x_sum_per_box[i];
        }
    }

    ParallelDescriptor::ReduceIntSum (&bin_count[0], pdf_num_bins);
    ParallelDescriptor::ReduceRealSum(&bin_x_sum[0], pdf_num_bins);

    // Normalization
    Real sum_ni_dxi = 0.0;
    for (unsigned int i = 0; i < pdf_num_bins; ++i) {
        Real dx_i = bin_edges[i+1] - bin_edges[i];
        sum_ni_dxi += bin_count[i] * dx_i;
    }

    // Normalization constant for N_i -> P_i
    Real pnorm = 1.0 / sum_ni_dxi;

    std::vector<Real> bin_x(pdf_num_bins);
    std::vector<Real> bin_p(pdf_num_bins);

    for (unsigned int i = 0; i < pdf_num_bins; ++i) {
        if (bin_count[i] > 0) {
            bin_x[i] = bin_x_sum[i] / Real(bin_count[i]);
            bin_p[i] = pnorm * Real(bin_count[i]);
        }
    }

    if (ParallelDescriptor::IOProcessor()) {

        std::ofstream pdf_output_file;
        std::stringstream filename;
        filename << field_name << "_pdf_nstep_" << std::setfill('0') << std::setw(5) << nStep;
        pdf_output_file.open(filename.str().c_str());
        pdf_output_file << std::scientific;
        for (unsigned int i = 0; i < pdf_num_bins; ++i) {
            pdf_output_file << std::setw(40) << bin_edges[i] << std::setw(15) << bin_count[i] << std::setw(15) << bin_x[i] << std::setw(15) << bin_p[i] << std::endl;
        }
        pdf_output_file.close();
    }

    return;
}
