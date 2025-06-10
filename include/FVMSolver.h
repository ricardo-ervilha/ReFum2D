#ifndef FVMSOLVER_H
#define FVMSOLVER_H

#include "Mesh.h"
#include "BoundaryCondition.h"

class FVMSolver {
    private:
        Mesh* mesh;
        double** A; //ponteiro
        double* _A; //ponteiro auxiliar
        double* b;
        double* skew;
        vector<BoundaryCondition*> boundaries;
        double gamma;

    public:
        FVMSolver(Mesh* mesh, BoundaryCondition *down, BoundaryCondition *right, BoundaryCondition *top, BoundaryCondition *left, double gamma);
        ~FVMSolver();
        void applyBoundariesConditions();
        void computeA();
        void computeb();
        void printA();
        void printB();
};

#endif