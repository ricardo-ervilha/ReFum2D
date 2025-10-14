#ifndef NSSOLVER_H
#define NSSOLVER_H

#include "pch.h"
#include <petscksp.h>

class Mesh;
class Cell;

class NSSolver{
    private:
        Mesh *mesh;
        
        Mat A_mom;
        Vec b_mom_x;
        Vec b_mom_y;

        Mat A_pc;
        Vec b_pc; 

        float mu;
        float rho;

        Vec u_face;
        Vec v_face;
        Vec p_face;

        Vec uc;
        Vec vc;
        Vec pc;

        Vec uc_aux;
        Vec vc_aux;
        Vec pc_aux;
        
        Vec mdotf;

        Vec ap;
        Vec mdotfcorr;
        Vec pcorr;
        Vec pfcorr;

        Vec ucorr;
        Vec vcorr;

        Vec wf;
        KSP ksp;
        void compute_wf();
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