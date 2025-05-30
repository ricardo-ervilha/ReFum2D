#include "../include/FVMSolver.h"

FVMSolver::FVMSolver(Mesh* mesh, BoundaryCondition *down, BoundaryCondition* right, BoundaryCondition* top, BoundaryCondition *left){
    this->mesh = mesh;
    this->boundaries.push_back(down);
    this->boundaries.push_back(right);
    this->boundaries.push_back(top);
    this->boundaries.push_back(left);

    int N = this->mesh->getNumCells();
    A = new double*[N];
    b = new double[N];

    for(int i = 0; i < N; i++){
        A[i] = new double[N];
        A[i] = 0;
        b[i] = 0;
    }
}

FVMSolver::~FVMSolver(){

}

void FVMSolver::applyBoundariesConditions(){
}

void FVMSolver::computeA(){

}

void FVMSolver::computeb(){

}