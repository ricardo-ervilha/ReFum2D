#include "FVMSolver.h"

int main(void){

    FVMSolver* solver = new FVMSolver();
    solver->readMesh("../inputs/malhaSimples.msh");

    return 0;
}