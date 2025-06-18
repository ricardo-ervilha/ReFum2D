#ifndef FVMSOLVER_H
#define FVMSOLVER_H

#include "Mesh.h"
#include "BoundaryCondition.h"
#include <iomanip>
#include <fstream>
#include <filesystem>
#include "../eigen-3.4.0/Eigen/Dense"

class FVMSolver {
    private:
        // Informações do problema: malha, constantes, etc.
        Mesh* mesh; // Estruturas de dados pré-processadas pela leitura da malha
        double (*gammafunc)(double, double); // Constante de difusão
        double (*sourcefunc)(double, double); // Termo fonte
        vector<BoundaryCondition*> boundaries;

        // Estruturas de dados relacionadas a montagem e solução do problema.
        Eigen::MatrixXd A;
        Eigen::VectorXd b;
        Eigen::VectorXd u;
        Eigen::VectorXd source; // pré-processa no FVMSolver os termos fontes das células.
        Eigen::VectorXd gammaf; // pré-processa no FVMSolver os gammas das faces;
        Eigen::VectorXd skew; // relacionado a angulosidade da malha, computa a difusão cruzada.

        void pre_processing(); // Função para realizar o pré-processamento de gamma, source.
        
        double phiv(Node* n); // Função auxiliar para calcular o valor de phi nos nós.

        double max_norm_diff(vector<double>& u1, vector<double>& u2);
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
        void apply_boundaries();

        /*Resolve sem considerar a difusão cruzada...*/
        void iterative_solver(double tol);

        // recebe o caminho e nome do arquivo
        void save_solution(string filepath);

        void summary() {this->mesh->mesh_summary();};

        double compute_error(double (*exact)(double, double));
        bool is_simetric();
};

#endif