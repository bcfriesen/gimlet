#ifndef DOGIMLETANALYSIS_H
#define DOGIMLETANALYSIS_H

void
do_analysis(const Real     omega_b,
            const Real     omega_m,
            const Real     omega_l,
            const Real     h,
            const Real     comoving_a,
            const MultiFab &state_data,
            const MultiFab &eos_data,
            const MultiFab &dm_density,
            const MultiFab &xmom,
            const MultiFab &ymom,
            const MultiFab &zmom,
            const Geometry &geom,
            const int      nStep);


#endif /* DOGIMLETANALYSIS_H */
