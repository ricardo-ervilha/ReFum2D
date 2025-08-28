#ifndef FVMSOLVER_H
#define FVMSOLVER_H

#include "pch.h"
#include "Mesh.h"
#include "BoundaryCondition.h"
#include <armadillo>

class FVMSolver
{
private:
    Mesh *mesh; // Estruturas de dados pré-processadas pela leitura da malha.

    double (*gammafunc)(double, double);           // Constante de difusão
    double (*sourcefunc)(double, double);          // Termo fonte
    double (*rhofunc)(double, double);             // Densidade do meio
    pair<double, double> (*ufunc)(double, double); // Vetor velocidade

    map<BoundaryLocation, BoundaryCondition *> boundaries; // Condições de contorno informados pelo usuário.

    // Estruturas de dados relacionadas a montagem e solução do problema.
    arma::sp_mat A;
    arma::vec b;
    arma::vec u_old; // u_n
    arma::vec u_new; // u_{n+1}
    arma::vec source;                       // pré-processa no FVMSolver os termos fontes das células.
    arma::vec gammaf;                       // pré-processa no FVMSolver os gammas das faces;
    arma::vec rhof;                         // pré-processa no FVMSolver os rhos das faces;
    arma::mat uf;                           // pré-processa no FVMSolver os us das faces;
    arma::vec b_corrected;                  // * Adiciona os termos vindos da iteração explícita.
    vector<pair<double, double>> gradients; // * vetor com os gradientes reconstruídos.

    void pre_processing(); // Função para realizar o pré-processamento dos dados do problema e colocar em ED's.

public:
    FVMSolver(Mesh *mesh, BoundaryCondition *bc1, BoundaryCondition *bc2, BoundaryCondition *bc3, BoundaryCondition *bc4, double (*g)(double, double), double (*rho)(double, double), pair<double, double> (*U)(double, double), double (*sourceTerm)(double, double));
    ~FVMSolver();

    void print_A()
    {
        cout << "#==================================================================================#" << endl;
        cout << A << endl;
        cout << "#==================================================================================#" << endl;
    };
    void print_b()
    {
        cout << "#==================================================================================#" << endl;
        cout << b << endl;
    };

    /* A partir da função definida pelo usuário inicializa a condição inicial do problema*/
    void set_initial_condition(double (*icFunc)(double, double));
    void transient(double dt); // * Calcula o transiente de fato.

    /* funções para calculo da convecção com upwind */
    void convection();
    void convection_of_cell(Cell *c);
    void correct_linear_upwind();

    /* funções para calculo da difusão sem considerar a não ortogonalidade */
    void diffusion_of_cell(Cell *c);
    void diffusion();

    /* função para reconstruir o gradiente*/
    void compute_gradients();

    /* função que utiliza a compute_gradients para corrigir a difusão cruzada */
    void compute_cross_diffusion();

    /* Resolve o sistema transiente salvando os snapshots em path */
    void TransientSolver(double t0, double tf, int snapshots, string path);

    /* Resolve o sistema estacionário */
    void SteadySolver(double tolerance);

    /* Usada quando se utiliza o SteadySolver */
    void save_solution(string filepath);

    void compute_error(double (*exact)(double, double));
};

#endif