#include "Mesh.h"
#include "NSSolver.h"

int main(int argc, char* argv[]) {
    Mesh m;
    m.read_mesh("../inputs/q30x30.msh");

    NSSolver solver(&m, 0.01, 1, 0, 0);

    // calcula coeficientes de difusão e convecção.
    solver.calculate_Df();
    solver.calculate_Gf();

    // calcula as equações de momento a partir do chute de pressão.
    solver.calculate_momentum();

    solver.export_solution("../outputs/u_star.vtk", 'u');
}
