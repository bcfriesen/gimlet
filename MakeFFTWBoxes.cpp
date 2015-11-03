#include <Geometry.H>
#include <BLFort.H>

/*
 * This routine generates a domain decomposition (characterized by a BoxList)
 * which FFTW can use to do distributed-memory DFTs. When doing
 * distributed-memory DFTs on multidimensional arrays, as we do here, FFTW
 * requires that the domain be striped along the last index in Fortran array
 * index notation, or along the first index in C-style index notation (the two
 * are functionally equivalent). Since FABs are laid out in column-major
 * (Fortran) order, we stripe along the last index (along the z-dimension).
 *
 * Although FFTW supports only one style of data striping (i.e., it does not
 * allow arbitrary box shapes and distributions like we have in BoxLib), it
 * nevertheless allows you to customize some of the finer points of how the
 * data are striped across processes. However you can also let it decide for
 * itself what the "best" data distribution is, which is what we do here.
 *
 * One final constraint is that we can have at most ONE box per process when we
 * do the DFT. We can't use a MFIter to iterate over bits and pieces of the
 * grid, doing partial DFTs on each Box; we must do the DFT on the entire grid
 * all at once. This has two consequences:
 *
 * 1.) Because FFTW stripes along the z-dimension, it is possible to have fewer
 * slabs than processes, in which cause the extra processes cannot participate
 * in the DFT. For example, if you run a simulation on a 1024^3 grid with 2048
 * processes, FFTW will generate at most 1024 slabs (each one cell thick) which
 * will be distributed on 1024 processes; the remaining 1024 processes will be
 * idle during the DFT execution call.
 *
 * 2.) Because we can't iterate over multiple Boxes per process during the DFT,
 * we cannot implement loop tiling with OpenMP threads in the same way we do
 * elsewhere in BoxLib. Happily, FFTW provides an internal threading model that
 * we can use instead (see the "Multi-threaded FFTW" section in the docs).
 */

BL_FORT_PROC_DECL(GET_FFTW_BOX_SIZES, get_fftw_box_sizes) (
        const intptr_t* domain_size,
        const int* comm,
        const intptr_t* local_z,
        const intptr_t* local_k_offset,
        const intptr_t* alloc_local
        );

void MakeFFTWBoxes(const Geometry& geom, BoxList& bl_fft, intptr_t& alloc_local) {

    const Box& domain = geom.Domain();
    const IntVect& domain_size_IV = domain.size();
    const int* const domain_size_int = domain_size_IV.getVect();
    const int comm = MPI_Comm_c2f(ParallelDescriptor::Communicator());

    // Most arguments to the Fortran interface to FFTW use this weird type
    // instead of regular int.
    intptr_t local_z, local_k_offset;
    Array<intptr_t> domain_size(3);
    for (unsigned int i = 0; i < 3; ++i) domain_size[i] = domain_size_int[i];

    BL_FORT_PROC_CALL(GET_FFTW_BOX_SIZES, get_fftw_box_sizes) (
            &domain_size[0],
            &comm,
            &local_z,
            &local_k_offset,
            &alloc_local);

    // Gather everybody's Box distribution so we can make the BoxArray
    const int nprocs = ParallelDescriptor::NProcs();
    Array<int> local_z_all(nprocs);
    Array<int> local_k_offset_all(nprocs);
    // MPI doesn't have a data type for intptr_t, so we cast to int. No idea if
    // this is safe on all architectures, but it seems to work on x86.
    const int local_z_int = int(local_z);
    const int local_k_offset_int = int(local_k_offset);
    MPI_Allgather(&local_z_int, 1, MPI_INT, &local_z_all[0], 1, MPI_INT, ParallelDescriptor::Communicator());
    MPI_Allgather(&local_k_offset_int, 1, MPI_INT, &local_k_offset_all[0], 1, MPI_INT, ParallelDescriptor::Communicator());

    // Construct Boxes as slabs. If we have more processes than total grid
    // points along the z-dimension, some processes will not get a slab and
    // will be idle during the call to FFTW.
    const IntVect& domain_smallEnd = domain.smallEnd();
    const IntVect& domain_bigEnd = domain.bigEnd();
    for (unsigned int i = 0; i < nprocs; ++i) {
        if (local_z_all[i] <= 0) continue; // If we have more procs than slabs, some of the slabs will have zero thickness.
        Box bx;
        bx.setSmall(0, domain_smallEnd[0]);
        bx.setBig  (0, domain_bigEnd  [0]);
        bx.setSmall(1, domain_smallEnd[1]);
        bx.setBig  (1, domain_bigEnd  [1]);
        bx.setSmall(2, domain_smallEnd[2] + local_k_offset_all[i]);
        bx.setBig  (2, domain_smallEnd[2] + local_k_offset_all[i] + local_z_all[i]-1);
        bl_fft.push_back(bx);
    }

    if (bl_fft.size() < nprocs) {
      if (ParallelDescriptor::IOProcessor()) {
        std::cerr << "WARNING: FFTW requested " << bl_fft.size() << " slabs, but " << nprocs << " processes are active. ";
        std::cerr << nprocs-bl_fft.size() << " processes will be idle during the DFT." << std::endl;
      }
    }
}
