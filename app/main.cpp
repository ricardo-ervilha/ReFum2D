#include "Mesh.h"
#include "NSSolver2.h"

int main(int argc, char* argv[]) {
    Mesh m;
    m.read_mesh("../inputs/q5x5.msh");

    NSSolver2 solver(&m, 1e-2, 1.0, 0, 0);
    
    for(int i = 0; i < 4; i++){
        solver.mom_links_and_sources();
        cout << "x-mom\n";
        solver.solve_x_mom(100000, 0.1, 1e-6);
        cout << "y-mom\n";
        solver.solve_y_mom(100000, 0.1, 1e-6);
        solver.face_velocity();
        cout << "pressure\n";
        solver.solve_pp(100, 1e-6);
        solver.uv_correct(0.7);
        solver.pres_correct(0.2);
    }

    
    solver.export_solution("../outputs/solution.vtk");
}
