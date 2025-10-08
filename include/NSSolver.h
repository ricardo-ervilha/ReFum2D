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
        std::map<BoundaryLocation, BoundaryCondition*> u_bcs;
        std::map<BoundaryLocation, BoundaryCondition*> v_bcs;
        std::map<BoundaryLocation, BoundaryCondition*> p_bcs;
        
        arma::sp_mat A_mom;
        arma::vec b_mom_x;
        arma::vec b_mom_y;

        arma::sp_mat A_pc;
        arma::vec b_pc; 

        float mu;
        float rho;

        vector<pair<BoundaryType, double>> u_face;
        vector<pair<BoundaryType, double>> v_face;
        vector<pair<BoundaryType, double>> p_face_bc;

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
        void compute_bc();
    public:
        NSSolver(Mesh *mesh, float mu, float rho, BoundaryCondition *bc1, BoundaryCondition* bc2, BoundaryCondition* bc3, BoundaryCondition *bc4, BoundaryCondition *bc5, BoundaryCondition* bc6, BoundaryCondition* bc7, BoundaryCondition *bc8, BoundaryCondition *bc9, BoundaryCondition* bc10, BoundaryCondition* bc11, BoundaryCondition *bc12);
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