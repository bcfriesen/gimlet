#include <iomanip>
#include <fstream>
#include <sstream>

#include <MultiFab.H>
#include <Geometry.H>
#include <BLFort.H>

BL_FORT_PROC_DECL(CALC_PDF2D, calc_pdf2d) (
        const Real* temperature,
        const Real* density,
        const int*  lo,
        const int*  hi,
        const int*  num_ghosts,
        const int*  tpdf_num_bins,
        const int*  dpdf_num_bins,
        const Real* t_bin_edges,
        const Real* d_bin_edges,
        const int*  bin_count,
        const Real* bin_t_sum,
        const Real* bin_d_sum);

void temperature_density_pdf2d (const MultiFab& temperature,
                                const MultiFab& density,
                                const Geometry &geom,
                                const Real mean_density,
                                const unsigned int nStep)
{

    MultiFab overdensity(density.boxArray(), 1, 0);
    overdensity.copy(density);
    overdensity.mult(1.0/mean_density);
    overdensity.plus(-1.0, 0, 1);

#warning TODO: make temperature PDF knobs run-time parameters
    const Real log10_t_min = 3.0;
    const Real log10_t_max = 8.0;
    const int tpdf_num_bins = 400;

#warning TODO: make density PDF knobs run-time parameters
    const Real log10_d_min = -2.0;
    const Real log10_d_max =  5.0;
    const int dpdf_num_bins = 400;

    const int temperature_num_ghosts = temperature.nGrow();

    // These are counts per *process*, not per *box*. We'll AllReduce() them
    // after iterating over all boxes.
    std::vector<int>  bin_count;
    std::vector<Real> bin_t_sum;
    std::vector<Real> bin_d_sum;
    std::vector<Real> bin_edges_t;
    std::vector<Real> bin_edges_d;

    const unsigned int num_bins = tpdf_num_bins*dpdf_num_bins;
    bin_count.resize(num_bins);
    bin_t_sum.resize(num_bins);
    bin_d_sum.resize(num_bins);
    for (unsigned int i = 0; i < num_bins; ++i) {
        bin_count[i] = 0;
        bin_t_sum[i] = 0.0;
        bin_d_sum[i] = 0.0;
    }

    bin_edges_t.resize(tpdf_num_bins+1);
    Real bin_width = (log10_t_max - log10_t_min) / Real(tpdf_num_bins);
    for (unsigned int i = 0; i <= tpdf_num_bins; ++i) {
        bin_edges_t[i] = std::pow(10.0, log10_t_min + Real(i)*bin_width);
    }

    bin_edges_d.resize(dpdf_num_bins+1);
    bin_width = (log10_d_max - log10_d_min) / Real(dpdf_num_bins);
    for (unsigned int i = 0; i <= dpdf_num_bins; ++i) {
        bin_edges_d[i] = std::pow(10.0, log10_d_min + Real(i)*bin_width);
    }

    for ( MFIter mfi(temperature); mfi.isValid(); ++mfi ) {
        const Box& bx = mfi.validbox();

        // These are counts per *box*. We'll glob them all together after
        // iterating over all boxes.
        std::vector<int>  bin_count_per_box(num_bins);
        std::vector<Real> bin_t_sum_per_box(num_bins);
        std::vector<Real> bin_d_sum_per_box(num_bins);

        BL_FORT_PROC_CALL(CALC_PDF2D, calc_pdf2d) (
                (temperature)[mfi].dataPtr(),
                (overdensity)[mfi].dataPtr(),
                bx.loVect(),
                bx.hiVect(),
                &temperature_num_ghosts,
                &tpdf_num_bins,
                &dpdf_num_bins,
                &bin_edges_t[0],
                &bin_edges_d[0],
                &bin_count_per_box[0],
                &bin_t_sum_per_box[0],
                &bin_d_sum_per_box[0]);

        // Accumulate per-box results into per-process arrays
        for (unsigned int i = 0; i < num_bins; ++i) {
             bin_count[i] += bin_count_per_box[i];
             bin_t_sum[i] += bin_t_sum_per_box[i];
             bin_d_sum[i] += bin_d_sum_per_box[i];
        }
    }

    ParallelDescriptor::ReduceIntSum (&bin_count[0], num_bins);
    ParallelDescriptor::ReduceRealSum(&bin_t_sum[0], num_bins);
    ParallelDescriptor::ReduceRealSum(&bin_d_sum[0], num_bins);

    // Normalization
    Real sum_ni_dai = 0.0;
    for (unsigned int i = 0; i < tpdf_num_bins; ++i) {
      for (unsigned int j = 0; j < dpdf_num_bins; ++j) {
        const unsigned idx = j + (i*dpdf_num_bins);
        Real dx_i = bin_edges_t[i+1] - bin_edges_t[i];
        Real dx_j = bin_edges_d[j+1] - bin_edges_d[j];
        sum_ni_dai += bin_count[idx] * dx_i * dx_j;
      }
    }

    // Normalization constant for N_i -> P_i
    Real pnorm = 1.0 / sum_ni_dai;

    std::vector<Real> bin_x(num_bins);
    std::vector<Real> bin_y(num_bins);
    std::vector<Real> bin_p(num_bins);

    for (unsigned int i = 0; i < num_bins; ++i) {
        if (bin_count[i] > 0) {
            bin_x[i] = bin_t_sum[i] / Real(bin_count[i]);
            bin_y[i] = bin_d_sum[i] / Real(bin_count[i]);
            bin_p[i] = pnorm * Real(bin_count[i]);
        }
    }

    if (ParallelDescriptor::IOProcessor()) {

        std::ofstream temperature_pdf_output_file;
        std::stringstream filename;
        filename << "rhot_pdf2d_nstep_" << std::setfill('0') << std::setw(5) << nStep;
        temperature_pdf_output_file.open(filename.str().c_str());
        temperature_pdf_output_file << std::scientific;
        for (unsigned int i = 0; i < dpdf_num_bins; ++i) {
          for (unsigned int j = 0; j < tpdf_num_bins; ++j) {
              const unsigned int idx = j + (i*tpdf_num_bins);
            temperature_pdf_output_file << std::setw(35) << bin_edges_d[i] << std::setw(35) << bin_edges_t[j] << std::setw(15) << bin_count[idx] << std::setw(15) << bin_x[idx] << std::setw(15) << bin_y[idx] << std::setw(15) << bin_p[idx] << std::endl;
          }
        }
        temperature_pdf_output_file.close();
    }

    return;
}
