#ifndef FVMSOLVER_H
#define FVMSOLVER_H

#include "pch.h"
#include "Mesh.h"
#include <armadillo>

class FVMSolver {
    private:
        // Informações do problema: malha, constantes, etc.
        Mesh* mesh; // Estruturas de dados pré-processadas pela leitura da malha
        double (*gammafunc)(double, double); // Constante de difusão
        double (*sourcefunc)(double, double); // Termo fonte
        vector<BoundaryCondition*> boundaries;

        // Estruturas de dados relacionadas a montagem e solução do problema.
        arma::mat A;
        arma::vec b;
        arma::vec u;
        arma::vec source; // pré-processa no FVMSolver os termos fontes das células.
        arma::vec gammaf; // pré-processa no FVMSolver os gammas das faces;
        arma::vec b_with_cd; // relacionado a angulosidade da malha, computa a difusão cruzada.
        vector<pair<double,double>> gradients; // vetor com os gradientes daquela iteração.
        
        void pre_processing(); // Função para realizar o pré-processamento de gamma, source.

        void print_matrix(vector<vector<double>>* m);
        void print_vector(vector<double> *v);
    
    public:    
        FVMSolver(string filepath, BoundaryCondition *down, BoundaryCondition *right, BoundaryCondition *top, BoundaryCondition *left, double (*g)(double, double), double (*sourceTerm)(double, double));
        ~FVMSolver();
        
        void print_A() { cout << A << endl;};
        void print_b() { cout << b << endl;};
        void print_Q() { cout << source << endl;};

        void assembly_A();
        void assembly_b();
        void compute_gradients();
        void compute_cross_diffusion();
        void apply_boundaries();

        /*Resolve sem considerar a difusão cruzada...*/
        void solve_system();

        // recebe o caminho e nome do arquivo
        void save_solution(string filepath);

        void summary() {this->mesh->mesh_summary();};

        void compute_error(double (*exact)(double, double));
        bool is_simetric();
};

#endif