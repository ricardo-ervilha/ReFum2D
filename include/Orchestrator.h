#ifndef ORCHESTRATOR_H
#define ORCHESTRATOR_H

#include "pch.h"
#include "BoundaryCondition.h"
class Mesh;
class Cell;

enum SolverType{
    Transient,
    Steady
};

class Orchestrator{
    private:
        Mesh* mesh;
        
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

        void readYamlAndRecoverVariables(string yaml_filepath);
    public:
        Orchestrator(string yaml_filepath);
        ~Orchestrator();
};

#endif