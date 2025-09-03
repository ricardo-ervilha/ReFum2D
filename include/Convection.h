#ifndef CONVECTION_H
#define CONVECTION_H

#include "pch.h"
#include <armadillo>

class Cell;
class FVMSolver;

class Convection{
    private:
        // & Função que recebe x, y para pegar o valor da densidade naquele vértice da malha.
        double (*rhofunc)(double, double); 

        // & Função que recebe x, y e pega o valor da velocidade no vértice.
        pair<double, double> (*Ufunc)(double, double);

        // # guarda o solver para recuperar malha e ed's
        FVMSolver* solver;

        void interpolate_properties();
        void control_volume_convection(vector<Cell*> cells, Cell* P); // calcula a convecção para cada volume de controle
    public:
        arma::vec rhof;
        arma::mat Uf;

        Convection(FVMSolver *s, double (*g)(double, double), pair<double, double> (*ufunc)(double, double));
        ~Convection();

        void assembleCoefficients();

        void linear_upwind();
};

#endif