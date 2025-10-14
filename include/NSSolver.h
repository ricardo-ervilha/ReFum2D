#ifndef NSSOLVER_H
#define NSSOLVER_H

#include "pch.h"
#include <armadillo>
#include "BoundaryCondition.h"

class Mesh;
class Cell;

class NSSolver{
    private:
        Mesh *mesh;

        vector<BoundaryCondition> bcsu;
        vector<BoundaryCondition> bcsv;
        vector<BoundaryCondition> bcsp;
        
        arma::sp_mat A_mom;
        arma::vec b_mom_x;
        arma::vec b_mom_y;

        arma::sp_mat A_pc;
        arma::vec b_pc; 

        float mu;
        float rho;

        vector<pair<BoundaryType,double>> u_boundary; 
        vector<pair<BoundaryType,double>> v_boundary; 
        vector<pair<BoundaryType,double>> p_boundary; 

        arma::vec u_face;
        arma::vec v_face;
        arma::vec p_face;

        arma::vec uc;
        arma::vec vc;
        arma::vec pc;

        arma::vec uc_aux;
        arma::vec vc_aux;
        arma::vec pc_aux;
        
        arma::vec mdotf;

        arma::vec ap;
        arma::vec mdotfcorr;
        arma::vec pcorr;
        arma::vec pfcorr;

        arma::vec ucorr;
        arma::vec vcorr;

        arma::vec wf;
        void compute_wf();
    public:
        NSSolver(Mesh *mesh, float mu, float rho, vector<BoundaryCondition> bcsU, vector<BoundaryCondition> bcsV, vector<BoundaryCondition> bcsP);
        ~NSSolver();

        void mom_links_and_sources(double lambda_v);
        void solve_x_mom();
        void solve_y_mom();
        void face_velocity();
        void solve_pp();
        void uv_correct();
        void pres_correct(double lambda_p);

        void compute_bcs();

        void export_solution(string filename);
};

#endif