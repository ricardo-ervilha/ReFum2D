#include "ConvectionBenchmark.h"
#include "FVMSolver.h"
#include "Mesh.h"
#include "Diffusion.h"
#include "Convection.h"
#include "Source.h"
#include "GradientReconstruction.h"

int main() {
    // Objetos automáticos
    Mesh m;
    m.read_mesh("../inputs/100x100.msh");

    BoundaryCondition downBC(DIRICHLET, DOWN, down);
    BoundaryCondition rightBC(DIRICHLET, RIGHT, right);
    BoundaryCondition topBC(DIRICHLET, TOP, top);
    BoundaryCondition leftBC(DIRICHLET, LEFT, left);

    FVMSolver solver(&m, &downBC, &rightBC, &topBC, &leftBC);

    Convection c(&solver, rho, velocity);
    GradientReconstruction g(&solver);

    c.assembleCoefficients();

    // for (int i = 0; i < 30; i++) {
        // reseta a contribuição explícita em b_aux
        solver.resetExplicit();
        
        // obtém os gradientes
        // g.reconstruct_gradients();

        // computa a parcela explícita do linear upwind
        // c.linear_upwind();

        // adiciona a contribuição de b em b_aux
        solver.addExplicitContribution();

        // resolve o sistema
        solver.solveSystem();
    // }


    // solver.compute_error(exact);
    solver.export_solution("../outputs/solution.vtk");

    return 0;
}
