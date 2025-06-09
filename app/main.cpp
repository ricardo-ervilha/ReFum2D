#include "../include/Mesh.h"
#include "../include/FVMSolver.h"
#include "../include/BoundaryCondition.h"

int main(void){

    Mesh* mesh = new Mesh();
    mesh->readMesh("../inputs/complexMesh.msh");
    mesh->meshSummary();
    BoundaryCondition down = BoundaryCondition("Dirichlet", 0.0);
    BoundaryCondition right = BoundaryCondition("Dirichlet", 0.0);
    BoundaryCondition top = BoundaryCondition("Dirichlet", 80.0);
    BoundaryCondition left = BoundaryCondition("Dirichlet", 0.0);
    FVMSolver sol = FVMSolver(mesh, &down, &right, &top, &left, 50);
    sol.computeA();
    return 0;
}