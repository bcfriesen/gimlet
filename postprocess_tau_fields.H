#ifndef POSTPROCESS_TAU_FIELDS_H
#define POSTPROCESS_TAU_FIELDS_H

#include <MultiFab.H>
#include <Geometry.H>

void postprocess_tau_fields(MultiFab& tau, const Geometry& geom, const int dir, const int nStep);

#endif
