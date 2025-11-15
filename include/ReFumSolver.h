#ifndef NSSOLVER_H
#define NSSOLVER_H

#include "pch.h"
#include <armadillo>
#include "BoundaryCondition.h"
#include "Orchestrator.h"
#include "eigen3/Eigen/Sparse"

class Mesh;
class Cell;

class ReFumSolver
{
private:
    Mesh *mesh;

    vector<BoundaryCondition> bcsu;
    vector<BoundaryCondition> bcsv;
    vector<BoundaryCondition> bcsp;

    Eigen::SparseMatrix<double> A;
    Eigen::VectorXd b_mom_x;
    Eigen::VectorXd b_mom_y;

    Eigen::VectorXd b_pc;

    float mu;
    float rho;
    float dt;

    vector<Eigen::Triplet<double>> triplets;
    Eigen::VectorXd diags; 

    vector<pair<BoundaryType, double>> u_boundary;
    vector<pair<BoundaryType, double>> v_boundary;
    vector<pair<BoundaryType, double>> p_boundary;

    Eigen::VectorXd u_face;
    Eigen::VectorXd v_face;
    Eigen::VectorXd p_face;

    Eigen::VectorXd uc;
    Eigen::VectorXd vc;
    Eigen::VectorXd pc;

    // & quando houver timestep, eles guardam os valores de condição inicial.
    Eigen::VectorXd uc_old;
    Eigen::VectorXd vc_old;
    Eigen::VectorXd pc_old;

    Eigen::VectorXd uc_aux;
    Eigen::VectorXd vc_aux;
    Eigen::VectorXd pc_aux;

    Eigen::VectorXd mdotf;

    Eigen::VectorXd ap;
    Eigen::VectorXd mdotfcorr;
    Eigen::VectorXd pcorr;
    Eigen::VectorXd pfcorr;

    SolverType solver;

    Eigen::VectorXd ucorr;
    Eigen::VectorXd vcorr;

    Eigen::VectorXd wf;
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