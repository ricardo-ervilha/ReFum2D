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
    // problem definitions
    string name;
    double rho;
    double nu;
    double reynolds;
    vector<BoundaryCondition> bcsu;
    vector<BoundaryCondition> bcsv;
    vector<BoundaryCondition> bcsp;

    string mshfile;

    // simple
    double lambda_uv;
    int non_corrections;
    int iterations_mom;
    double tolerance_mom;

    double lambda_p;
    int iterations_pc;
    double tolerance_pc;

    double utol, vtol, ptol;

    bool save_iterations;
    string exportfolder;
    
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

    // strings
    const std::string& get_name() const { return name; }
    const std::string& get_mshfile() const { return mshfile; }
    const std::string& get_exportfolder() const { return exportfolder; }

    // propriedades físicas
    double get_rho() const { return rho; }
    double get_nu() const { return nu; }
    double get_reynolds() const { return reynolds; }

    // condições de contorno
    const std::vector<BoundaryCondition>& get_bcsu() const { return bcsu; }
    const std::vector<BoundaryCondition>& get_bcsv() const { return bcsv; }
    const std::vector<BoundaryCondition>& get_bcsp() const { return bcsp; }

    // SIMPLE
    double get_lambda_uv() const { return lambda_uv; }
    int get_non_corrections() const { return non_corrections; }
    int get_iterations_mom() const { return iterations_mom; }
    double get_tolerance_mom() const { return tolerance_mom; }

    double get_lambda_p() const { return lambda_p; }
    int get_iterations_pc() const { return iterations_pc; }
    double get_tolerance_pc() const { return tolerance_pc; }

    // tolerâncias
    double get_utol() const { return utol; }
    double get_vtol() const { return vtol; }
    double get_ptol() const { return ptol; }

    // controle
    bool get_save_iterations() const { return save_iterations; }

    std::vector<BoundaryCondition>& get_bcsu() {
         return bcsu;
    }

    std::vector<BoundaryCondition>& get_bcsv() {
        return bcsv;
    }

    std::vector<BoundaryCondition>& get_bcsp() {
        return bcsp;
    }

};

#endif