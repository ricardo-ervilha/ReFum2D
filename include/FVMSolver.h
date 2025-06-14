#ifndef FVMSOLVER_H
#define FVMSOLVER_H

#include "Mesh.h"
#include "BoundaryCondition.h"
#include <iomanip>
#include <fstream>
#include <filesystem>

class FVMSolver {
    private:
        Mesh* mesh;
        double** A; //ponteiro
        double* _A; //ponteiro auxiliar
        double* b;
        double* skew;
        vector<BoundaryCondition*> boundaries;
        double gamma; // constante de difusividade
        double *u; // vetor solução (acompanha indexação das células, ou seja, cada entrada sua é a solução correspondente do valor de u para o centroid daquele determinado volume de controle)

        double *Qs;
        
        double calculateQ(double x, double y);
        double phiv(Node* n); // função auxiliar para calcular o valor de phi nos nós~
        double* copyVector(double *M, int N); 
        double max_norm_difference(double* u, double* u_old, int N);
    public:
        FVMSolver(Mesh* mesh, BoundaryCondition *down, BoundaryCondition *right, BoundaryCondition *top, BoundaryCondition *left, double gamma);
        ~FVMSolver();
        void applyBoundariesConditions();
        void computeA();
        void computeb();
        void printA();
        void printB();
        void computeQs();
        void GaussSeidel(double tol);
        void saveSolution(string filename);
        void computeErrorExact();
};

#endif