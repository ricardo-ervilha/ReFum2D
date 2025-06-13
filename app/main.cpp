#include "../include/Mesh.h"
#include "../include/FVMSolver.h"
#include "../include/BoundaryCondition.h"

int main(void){

    Mesh* mesh = new Mesh();
    mesh->readMesh("../inputs/mesh1.msh");
    mesh->meshSummary();
    BoundaryCondition down = BoundaryCondition("Dirichlet", 0.0);
    BoundaryCondition right = BoundaryCondition("Dirichlet", 0.0);
    BoundaryCondition top = BoundaryCondition("Dirichlet", 0.0);
    BoundaryCondition left = BoundaryCondition("Dirichlet", 0.0);
    FVMSolver sol = FVMSolver(mesh, &down, &right, &top, &left, 1);
    double tol = 1e-3;
    sol.computeA();
    sol.computeb();
    // sol.applyBoundariesConditions();
    sol.GaussSeidel(tol);
    sol.saveSolution();
    sol.computeErrorExact();
    return 0;
}