#ifndef FVMSOLVER_H
#define FVMSOLVER_H

#include "pch.h"
#include "BoundaryCondition.h"
#include <armadillo>

class Mesh;

/**
 * Comentários:
 * TODO: Indicar tarefas pendentes ou algo que precisa ser implementado.
 * ! Marcar algo muito importante ou crítico no código.
 * + Explicações de contas, fórmulas ou raciocínio matemático dentro do código.
 * * Comentários gerais ou explicativos sobre o código, documentação rápida.
 */

class FVMSolver
{
public:
    // # Ponteiro para a malha, podendo acessar células, arestas, nós, etc.
    Mesh *mesh; 

    // # Condições de contorno informados pelo usuário.
    map<BoundaryLocation, BoundaryCondition *> boundaries; 

    // # Sistema linear.
    arma::sp_mat A; // & Matriz A, de A x = b
    arma::vec b;    // & Matriz b sem correções explícitas, de A x = b
    arma::vec b_aux; // & VETOR AUXILIAR QUE TERÁ AS CORREÇÕES
    arma::vec u_old; // & u^{n} ou valor da solução no passo de tempo anterior
    arma::vec u_new; // & u^{n+1} ou valor da solução no próximo passo de tempo 
    arma::mat gradients; // & Armazena os gradientes reconstruídos com LSQ
    
    FVMSolver(Mesh *mesh, BoundaryCondition *bc1, BoundaryCondition *bc2, BoundaryCondition *bc3, BoundaryCondition *bc4);
    ~FVMSolver();

    // ! Computa o erro da solução se houver uma exata passando
    void compute_error(double (*exact)(double, double), string suffix_filepath);

    // ! Exporta a solução em formato .vtk
    void export_solution(string filepath);

    // Funções para ajudar no loop da solução
    void resetExplicit() { b_aux.zeros(); }
    void addExplicitContribution() { b_aux += b; }
    void solveSystem() { u_new = arma::spsolve(A, b_aux); }

    void set_initial_condition(double (*icFunc)(double, double));
};

#endif