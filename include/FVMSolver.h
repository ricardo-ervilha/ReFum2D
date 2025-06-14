#ifndef FVMSOLVER_H
#define FVMSOLVER_H

#include "Mesh.h"
#include "BoundaryCondition.h"
#include <iomanip>
#include <fstream>
#include <filesystem>

class FVMSolver {
    private:
        // Informações do problema: malha, constantes, etc.
        Mesh* mesh; // Estruturas de dados pré-processadas pela leitura da malha
        double (*gammafunc)(double, double); // Constante de difusão
        double (*sourcefunc)(double, double); // Termo fonte
        vector<BoundaryCondition*> boundaries;

        // Estruturas de dados relacionadas a montagem e solução do problema.
        vector<vector<double>> A;
        vector<double> b;
        vector<double> u;
        vector<double> source; // pré-processa no FVMSolver os termos fontes das células.
        vector<double> gamma; // pré-processa no FVMSolver os gammas das células;
        vector<double> skew; // relacionado a angulosidade da malha, computa a difusão cruzada.

        void pre_processing(); // Função para realizar o pré-processamento de gamma, source.
        
        double phiv(Node* n); // Função auxiliar para calcular o valor de phi nos nós.

        double max_norm_diff(vector<double>& u1, vector<double>& u2);
    
    public:    
        FVMSolver(Mesh* mesh, BoundaryCondition *down, BoundaryCondition *right, BoundaryCondition *top, BoundaryCondition *left, double (*g)(double, double), double (*sourceTerm)(double, double));
        ~FVMSolver();
        
        void print_matrix(vector<vector<double>>* m);
        void print_vector(vector<double> *v);

        void assembly_A();
        void assembly_b();

        // recebe o caminho e nome do arquivo
        void save_solution(string filepath);
};

#endif