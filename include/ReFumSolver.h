#ifndef NSSOLVER_H
#define NSSOLVER_H

#include "pch.h"
#include <armadillo>
#include "BoundaryCondition.h"
#include "Orchestrator.h"

class Mesh;
class Cell;

class ReFumSolver
{
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
    float dt;

    vector<pair<BoundaryType, double>> u_boundary;
    vector<pair<BoundaryType, double>> v_boundary;
    vector<pair<BoundaryType, double>> p_boundary;

    arma::vec u_face;
    arma::vec v_face;
    arma::vec p_face;

    arma::vec uc;
    arma::vec vc;
    arma::vec pc;

    // & quando houver timestep, eles guardam os valores de condição inicial.
    arma::vec uc_old;
    arma::vec vc_old;
    arma::vec pc_old;

    arma::vec uc_aux;
    arma::vec vc_aux;
    arma::vec pc_aux;

    arma::vec mdotf;

    arma::vec ap;
    arma::vec mdotfcorr;
    arma::vec pcorr;
    arma::vec pfcorr;

    SolverType solver;

    arma::vec ucorr;
    arma::vec vcorr;

    arma::vec wf;
    void compute_wf();
    void compute_bcs_first();
    void compute_bcs_repeat();

    // exportar velocidade e pressão em .vtk
    void export_velocity(string filename);
    void export_pressure(string filename);

    vector<vector<int>> cnsigns;
    vector<double>      careas;
    vector<pair<double,double>> ccentroids;
    vector<double>      flengths;
    vector<double>      fdfs;
    vector<pair<double,double>> fmiddles;
    vector<pair<double,double>> fnormals;
    vector<pair<int,int>> flftcs;
    vector<bool> fboundaryfaces;
    vector<vector<int>> idFacesFromCell;

public:
    ReFumSolver(Mesh *mesh, float mu, float rho, vector<BoundaryCondition> bcsU, vector<BoundaryCondition> bcsV, vector<BoundaryCondition> bcsP, SolverType solver);
    ~ReFumSolver();

    void set_initial_condition(function<double(double, double)> u_func, function<double(double, double)> v_func, function<double(double, double)> p_func);

    // simple
    void mom_links_and_sources(double lambda_v);
    void solve_x_mom();
    void solve_y_mom();
    void face_velocity();
    void solve_pp(bool sing_matrix);
    void uv_correct();
    void pres_correct(double lambda_p);

    void calculate_exact_solution_and_compare();

    // solver
    void TRANSIENTE_SIMPLE(string problem, string filepath, int num_simple_iterations, double lambda_uv, double lambda_p, int n_steps, double tf, bool pressure_correction_flag);
    void STEADY_SIMPLE(string problem, string filepath, int num_simple_iterations, double lambda_uv, double lambda_p, bool pressure_correction_flag);
};

#endif