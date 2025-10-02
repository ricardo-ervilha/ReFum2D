#ifndef NSSOLVER2_H
#define NSSOLVER2_H

#include "pch.h"
#include <armadillo>

class Mesh;
class Cell;

class NSSolver2{
    private:
        Mesh *mesh;
        
        
        // Viscosidade dinâmica
        float mu;
        // Densidade
        float rho;
        // termos fonte
        float source_x, source_y;

        arma::vec u_face;
        arma::vec v_face;

        arma::vec uc;
        arma::vec vc;
        arma::vec pc;
        
        arma::vec p_face;
        arma::vec mdotf;

        arma::vec scx;
        arma::vec scy;
        arma::vec ap;
        arma::vec res;
        arma::mat anb; // n linhas x 4 colunas

        arma::vec sc_p;
        arma::vec ap_p;
        arma::vec res_p;
        arma::mat anb_p;

        arma::vec mdotfcorr;
        arma::vec pcorr;
        arma::vec pfcorr;

        arma::vec ucorr;
        arma::vec vcorr;
    public:
        NSSolver2(Mesh *mesh, float mu, float rho, float source_x, float source_y);
        ~NSSolver2();

        void mom_links_and_sources();
        void solve_x_mom(int iter_mom, double rin_uv, double tol_inner);
        void solve_y_mom(int iter_mom, double rin_uv, double tol_inner);
        void face_velocity();
        void solve_pp(int iter_pp, double tol_inner);
        void uv_correct(double relax_uv);
        void pres_correct(double relax_p);

        void export_solution(string filepath);

};

#endif