#ifndef TEMPERATURE_DENSITY_PDF2D_H
#define TEMPERATURE_DENSITY_PDF2D_H

void temperature_density_pdf2d (const MultiFab& temperature, const MultiFab& density, const Geometry &geom, const Real mean_density, const unsigned int nStep);

#endif
