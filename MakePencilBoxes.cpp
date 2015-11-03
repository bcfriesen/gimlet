#include <Box.H>
#include <Geometry.H>

BoxList MakePencilBoxes(const Geometry& geom, const int dir) {

    const Box& problem_domain = geom.Domain();
    const int grid_nx = problem_domain.length(0);
    const int grid_ny = problem_domain.length(1);
    const int grid_nz = problem_domain.length(2);

    BoxList box_list;

    // pencils point in x-direction
    if (dir == 0) {
        for (int k = 0; k < grid_nz; ++k) {
            for (int j = 0; j < grid_ny; ++j) {
                Box pencil;
                pencil.setSmall(0, 0);
                pencil.setBig(0, grid_nx-1);
                pencil.setSmall(1, j);
                pencil.setBig(1, j);
                pencil.setSmall(2, k);
                pencil.setBig(2, k);
                box_list.push_back(pencil);
            }
        }
    // pencils point in y-direction
    } else if (dir == 1) {
        for (int k = 0; k < grid_nz; ++k) {
            for (int i = 0; i < grid_nx; ++i) {
                Box pencil;
                pencil.setSmall(0, i);
                pencil.setBig(0, i);
                pencil.setSmall(1, 0);
                pencil.setBig(1, grid_ny-1);
                pencil.setSmall(2, k);
                pencil.setBig(2, k);
                box_list.push_back(pencil);
            }
        }
    // pencils point in z-direction
    } else if (dir == 2) {
        for (int j = 0; j < grid_ny; ++j) {
            for (int i = 0; i < grid_nx; ++i) {
                Box pencil;
                pencil.setSmall(0, i);
                pencil.setBig(0, i);
                pencil.setSmall(1, j);
                pencil.setBig(1, j);
                pencil.setSmall(2, 0);
                pencil.setBig(2, grid_nz-1);
                box_list.push_back(pencil);
            }
        }
    }

    return box_list;
}
