#include "WolframAlpha.h"
#include "FVMSolver.h"
#include "Mesh.h"
#include "Diffusion.h"
#include "Convection.h"
#include "Source.h"
#include "GradientReconstruction.h"

int main() {
    // Objetos automáticos
    Mesh m;
    m.read_mesh("../inputs/quadStretchRefined.msh");

    BoundaryCondition downBC(NEUMANN, DOWN, down);
    BoundaryCondition rightBC(DIRICHLET, RIGHT, right);
    BoundaryCondition topBC(NEUMANN, TOP, top);
    BoundaryCondition leftBC(DIRICHLET, LEFT, left);

    FVMSolver solver(&m, &downBC, &rightBC, &topBC, &leftBC);

    Diffusion d(&solver, gamma);
    Source s(&solver, source);
    Convection c(&solver, rho, velocity);
    GradientReconstruction g(&solver);

    d.assembleCoefficients();
    s.assemblyCoefficients();
    c.assembleCoefficients();

    for (int i = 0; i < 10; i++) {
        // reseta a contribuição explícita em b_aux
        solver.resetExplicit();
        
        // obtém os gradientes
        g.reconstruct_gradients();

        // corrige a difusão cruzada
        d.cross_diffusion();

        // computa a parcela explícita do linear upwind
        c.linear_upwind();

        // adiciona a contribuição de b em b_aux
        solver.addExplicitContribution();

        // resolve o sistema
        solver.solveSystem();
    }


    solver.compute_error(exact);
    solver.export_solution("../outputs/benchmarkConvection.vtk");

    return 0;
}
