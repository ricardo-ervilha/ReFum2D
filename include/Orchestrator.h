#ifndef ORCHESTRATOR_H
#define ORCHESTRATOR_H

#include "pch.h"
#include "BoundaryCondition.h"
#include <functional>
class Mesh;
class Cell;

enum SolverType
{
    Transient,
    Steady
};

class Orchestrator
{
private:
    string problem;
    string mesh_path;

    vector<BoundaryCondition> bcsu;
    vector<BoundaryCondition> bcsv;
    vector<BoundaryCondition> bcsp;

    double rho;
    double nu;

    SolverType st;
    int n_steps;
    int tf;

    double lambda_uv, lambda_p;
    int iterations;
    string export_path;

    function<double(double, double)> u_ic;
    function<double(double, double)> v_ic;
    function<double(double, double)> p_ic;

    map<string, function<double(double, double)>> name_to_function_object;

    /* Funções necessárias para valor e condição inicial serão definidas aqui. */
    static double one(double x, double y) { return 1.0; };
    static double zero(double x, double y) { return 0.0; };
    
    static double backward_facing_step(double x, double y)
    {
        if (y >= 0.5) // metade superior
            return -16 * pow((y - 0.75), 2) + 1;
        else // metade inferior.
            return 0.0;
    };
    static double flow_over_cylinder_benchmark_a(double x, double y)
    {
        return 4 * 0.3 * y * (0.41 - y) / (0.41 * 0.41); // Re = 20
    }
    static double flow_over_cylinder_benchmark_b(double x, double y){
        return 4 * 1.5 * y * (0.41 - y) / (0.41 * 0.41); // U_m multiplicou por 5 => Re = 100
    }
    static double kov_u(double x, double y){
        double lambda = -1.81009812;
        return 1 - exp(lambda*x) * cos(2*M_PI*y);
    }
    static double kov_v(double x, double y){
        double lambda = -1.81009812;
        return (lambda/(2*M_PI))*exp(lambda*x)*sin(2*M_PI*y);
    }
    static double kov_p(double x, double y){
        double lambda = -1.81009812;
        return 0.5 * (1 - exp(2 * lambda * x));
    }



public:
    Orchestrator();
    ~Orchestrator();
    void readYamlAndRecoverVariables(string yaml_filepath);

    std::string get_problem() const { return problem; }
    std::string get_mesh_path() const { return mesh_path; }

    std::vector<BoundaryCondition> get_ubcs() const { return bcsu; }
    std::vector<BoundaryCondition> get_vbcs() const { return bcsv; }
    std::vector<BoundaryCondition> get_pbcs() const { return bcsp; }

    double get_rho() const { return rho; }
    double get_nu() const { return nu; }

    SolverType get_solver_type() const { return st; }
    int get_n_steps() const { return n_steps; }
    int get_tf() const { return tf; }

    double get_lambda_uv() const { return lambda_uv; }
    double get_lambda_p() const { return lambda_p; }

    int get_iterations() const { return iterations; }
    std::string get_export_path() const { return export_path; }

    function<double(double, double)> get_uic() {return u_ic;}
    function<double(double, double)> get_vic() {return v_ic;}
    function<double(double, double)> get_pic() {return p_ic;}
};

#endif