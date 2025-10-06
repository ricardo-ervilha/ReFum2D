#ifndef NSSOLVER_H
#define NSSOLVER_H

#include "pch.h"
#include <armadillo>

class Mesh;
class Cell;

class NSSolver{
    private:
        Mesh *mesh;
        
        arma::mat A_mom;
        arma::vec b_mom_x;
        arma::vec b_mom_y;

        arma::mat A_pc;
        arma::vec b_pc; 

        float mu;
        float rho;

        arma::vec u_face;
        arma::vec v_face;
        arma::vec p_face;

        arma::vec uc;
        arma::vec vc;
        arma::vec pc;
        
        arma::vec mdotf;

        arma::vec ap;
        arma::vec mdotfcorr;
        arma::vec pcorr;
        arma::vec pfcorr;

        arma::vec ucorr;
        arma::vec vcorr;
    public:
        NSSolver(Mesh *mesh, float mu, float rho);
        ~NSSolver();

        void mom_links_and_sources(double lambda_v);
        void solve_x_mom();
        void solve_y_mom();
        void face_velocity();
        void solve_pp();
        void uv_correct();
        void pres_correct(double lambda_p);

        void export_solution(string filename);
};

#endif