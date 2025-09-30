#include "Mesh.h"
#include "NSSolver.h"

int main(int argc, char* argv[]) {
    Mesh m;
    m.read_mesh("../inputs/q5x5.msh");

    NSSolver solver(&m, 1, 1, 0, 0);

    solver.calculate_Df();
    solver.calculate_Gf();
    
    for(int i = 0; i < 10; i++){
        solver.calculate_momentum();
        solver.interpolate_momentum();
        solver.pressure_correction_poisson();
        solver.correct_variables(0.7, 0.2);
    }
    
    solver.export_solution("../outputs/u_APPROX.vtk", 'u');
    solver.export_solution("../outputs/v_APPROX.vtk", 'v');
    solver.export_solution("../outputs/p_APPROX.vtk", 'p');
    // exatos..
    // solver.export_solution("../outputs/u_exact.vtk", 'x');
    // solver.export_solution("../outputs/v_exact.vtk", 'y');
    // solver.export_solution("../outputs/p_exact.vtk", 'z');
}
