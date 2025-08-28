#ifndef DIFFUSION_H
#define DIFFUSION_H

class FVMSolver;
class Cell;
#include <armadillo>
#include "pch.h"

class Diffusion{
    // & Função que recebe x, y para pegar o valor do coeficiente de difusão naquele vértice da malha.
    double (*gammafunc)(double, double); 

    // # guarda o solver para recuperar malha e ed's
    FVMSolver* solver;

    // & Interpola os gammas em gammaf.
    void interpolate_gammas();

    // + Wrapper para calcular para uma célula só.
    void control_volume_non_orthogonal_diffusion(vector<Cell*> cells, Cell *P);

    // & Computa o vetor n1
    pair<double,double> compute_n1(pair<double,double>& p, pair<double,double>& n, pair<double,double>& normal);

    public:
        // # Irá armazenar o valor de gamma em cada face interpolado 
        arma::vec gammaf;

        Diffusion(FVMSolver* s, double (*g)(double, double));
        ~Diffusion();

        // ! Irá adicionar os coeficientes p/ cada VC sem considerar a questão de não ortogonalidade
        void assembleCoefficients();

        // ! Irá adicionados os coeficientes usando gradientes reconstruídos para consertar a não ortogonalidade da malha
        void cross_diffusion();
};

#endif