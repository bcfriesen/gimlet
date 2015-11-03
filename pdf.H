#ifndef PDF_H
#define PDF_H

void pdf (const MultiFab& mf,
          const Geometry& geom,
          const unsigned int nStep,
          const std::string field_name,
          const Real log10_min,
          const Real log10_max,
          const int pdf_num_bins);

#endif
