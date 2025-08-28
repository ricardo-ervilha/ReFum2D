#include "WolframAlpha.h"
#include "FVMSolver.h"
#include "Diffusion.h"
#include "Convection.h"
#include "Source.h"

int main(void){
    // Instancia objeto da malha
    Mesh* m = new Mesh();
    
    // Leitura da malha e pré-processamento.
    m->read_mesh("../inputs/quadStretchRefined.msh");

    // Condições de contorno
    BoundaryCondition* downBC = new BoundaryCondition(NEUMANN, DOWN, down);
    BoundaryCondition* rightBC = new BoundaryCondition(DIRICHLET, RIGHT, right);
    BoundaryCondition* topBC = new BoundaryCondition(NEUMANN, TOP, top);
    BoundaryCondition* leftBC = new BoundaryCondition(DIRICHLET, LEFT, left);
    FVMSolver* solver = new FVMSolver(m, downBC, rightBC, topBC, leftBC);

    Diffusion* d = new Diffusion(solver, gamma);
    Source* s = new Source(solver, source);
    Convection* c = new Convection(solver, rho, velocity);

    d->assembleCoefficients();
    s->assemblyCoefficients();
    c->assembleCoefficients();

    solver->SteadySolver(d, true, 10);

    solver->compute_error(exact);

    solver->export_solution("../outputs/benchmarkDiffusion.vtk");

    /* Limpando ponteiros antes declarados */
    delete m;
    delete downBC;
    delete rightBC;
    delete topBC;
    delete leftBC;
    delete solver;
    delete d;
    delete s;

    return 0;
}