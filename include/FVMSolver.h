#ifndef FVMSOLVER_H
#define FVMSOLVER_H

#include "pch.h"
#include "Mesh.h"
#include "BoundaryCondition.h"
#include <armadillo>

class FVMSolver {
    private:
        Mesh* mesh; // Estruturas de dados pré-processadas pela leitura da malha.
        
        double (*gammafunc)(double, double); // Constante de difusão
        double (*sourcefunc)(double, double); // Termo fonte
        double (*rhofunc)(double, double); // Densidade do meio
        pair<double,double> (*ufunc)(double, double); // Vetor velocidade 
        
        map<BoundaryLocation, BoundaryCondition*> boundaries; // Condições de contorno informados pelo usuário.

        // Estruturas de dados relacionadas a montagem e solução do problema.
        arma::sp_mat A;
        arma::vec b;
        arma::vec u;
        arma::vec source; // pré-processa no FVMSolver os termos fontes das células.
        arma::vec gammaf; // pré-processa no FVMSolver os gammas das faces;
        arma::vec rhof; // pré-processa no FVMSolver os rhos das faces;
        arma::mat uf;    // pré-processa no FVMSolver os us das faces; 
        arma::vec b_with_cd; // relacionado a angulosidade da malha, computa a difusão cruzada.
        vector<pair<double,double>> gradients; // vetor com os gradientes daquela iteração.
        
        void pre_processing(); // Função para realizar o pré-processamento dos dados do problema e colocar em ED's.
        
        void diffusion(Cell* cell, Edge* edge, int nsign); // aplica a difusão na A (diferenças centradas)
        void convection(Cell* cell, Edge* edge, int nsign); // aplica a convecção na A (Upwind)
    public:    
        FVMSolver(Mesh* mesh, BoundaryCondition *bc1, BoundaryCondition *bc2, BoundaryCondition *bc3, BoundaryCondition *bc4, double (*g)(double, double), double (*rho)(double,double), pair<double,double> (*U)(double, double), double (*sourceTerm)(double, double));
        ~FVMSolver();
        
        void print_A() { 
            cout << "#==================================================================================#" << endl;
            cout << A << endl;
            cout << "#==================================================================================#" << endl;
        };
        void print_b() { 
            cout << "#==================================================================================#" << endl;
            cout << b << endl;
        };

        void assembly_A();
        void assembly_b();
        void compute_gradients();
        void compute_cross_diffusion();

        void solve_system(double tolerance);

        void save_solution(string filepath);

        void compute_error(double (*exact)(double, double));
};

#endif