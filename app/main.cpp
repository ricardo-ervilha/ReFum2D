#include "../include/Mesh.h"
#include "../include/FVMSolver.h"
#include "../include/BoundaryCondition.h"

int main(void){
    //mesh 1 => 120 triangles (0.15) : erro máx: 0.13
    //mesh 2 => 944 triangles (0.05) : erro máx: 0.23
    //mesh 3 => 1478 triangles (0.04): erro máx: 0.17
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
    sol.saveSolution("../outputs/mesh1.vtk");
    sol.computeErrorExact();
    return 0;
}