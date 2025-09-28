#ifndef NSSOLVER_H
#define NSSOLVER_H

#include "pch.h"
#include <armadillo>

class Mesh;
class Cell;

class NSSolver{
    private:
        Mesh *mesh;
        
        // + Momento
        // * Matriz A do sistema A x = b.
        arma::mat A_mom;
        // * Vetor coluna b do sistema A x = b.
        arma::vec b_mom; 

        // + Correção da pressão
        // * Matriz A do sistema A x = b.
        arma::mat A_pc;
        // * Vetor coluna b do sistema A x = b. 
        arma::vec b_pc; 

        // * Velocidade horizontal (u) atualizada [u = u* + u'].
        arma::vec u;
        // * Velocidade horizontal (u*) aproximada.
        arma::vec u_star;

        // * Velocidade vertical (v) atualizada [v = v* + v'].
        arma::vec v;
        // * Velocidade vertical (v*) aproximada.
        arma::vec v_star;

        // * Pressão (p) atualizada [p = p* + p'].
        arma::vec p;
        // * Pressão (p*) aproximada.
        arma::vec p_star;
        // * Correção da pressão (p1).
        arma::vec p_prime;

        // * Velocidade horizontal (u) na face.
        arma::vec u_face;
        // * Velocidade vertical (v) na face.
        arma::vec v_face;

        // * Coeficiente de difusão de cada volume de controle.
        arma::vec Df;
        // * Coeficiente de convecção de cada volume de controle.
        arma::vec Gf;

        arma::vec u_face_prime;
        arma::vec v_face_prime;

        // Viscosidade dinâmica
        float mu;
        // Densidade
        float rho;
        // termos fonte
        float source_x, source_y;

        arma::vec a;
    
    public:
        NSSolver(Mesh *mesh, float mu, float rho, float source_x, float source_y);
        ~NSSolver();

        // Função para calcular o coeficiente de difusão Df. Constantes durante toda a execução.
        void calculate_Df();
        // Função para calcular o coeficiente de convecção Gf. Serão recalculados a cada iteração externa.
        void calculate_Gf();

        // Calcula o momento para u e v.
        void calculate_momentum();
        // Usado no cálculo do momento para calcular os coeficientes de A.
        void calculate_A_mom();
        // Usado no calculo do momento para calcular os coeficientes de b. Leva em consdeiração a variavel 
        void calculate_b_mom_x();
        void calculate_b_mom_y();

        // Interpola o momento 
        void interpolate_momentum();

        // reconstroi gradientes de pressão;
        pair<double,double> reconstruct_pressure_gradients(Cell *c);

        void pressure_correction_poisson();

        void correct_variables();

        void export_solution(string filepath);

        void print_Df() {cout << this->Df << endl;};
        void print_Gf() {cout << this->Gf << endl;};

};

#endif