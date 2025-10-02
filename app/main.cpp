#include "Mesh.h"
#include "NSSolver2.h"

int main(int argc, char* argv[]) {
    Mesh m;
    m.read_mesh("../inputs/lid_driven_60.msh");

    NSSolver2 solver(&m, 1e-3, 1e3, 0, 0);
    
    for(int i = 0; i < 10; i++){
        solver.mom_links_and_sources();
        solver.solve_x_mom(100000, 0.1, 1e-6);
        solver.solve_y_mom(100000, 0.1, 1e-6);
        solver.face_velocity();
        solver.solve_pp(100, 1e-6);
        solver.uv_correct(0.8);
        solver.pres_correct(0.15);
    }

    
    solver.export_solution("../outputs/solution.vtk");
}
