#ifndef FVMSOLVER_H
#define FVMSOLVER_H

#include "Mesh.h"
#include "BoundaryCondition.h"

class FVMSolver {
    private:
        Mesh* mesh;
        double** A;
        double* b;
        vector<BoundaryCondition*> boundaries;

    public:
        FVMSolver(Mesh* mesh, BoundaryCondition *down, BoundaryCondition *right, BoundaryCondition *top, BoundaryCondition *left);
        ~FVMSolver();
        void applyBoundariesConditions();
        void computeA();
        void computeb();
};

#endif