#include "Orchestrator.h"
#include <yaml-cpp/yaml.h>

Orchestrator::Orchestrator(){
    name_to_function_object.emplace("one", one);
    name_to_function_object.emplace("zero", zero);
    name_to_function_object.emplace("bfs", backward_facing_step);
    name_to_function_object.emplace("foc_a", flow_over_cylinder_benchmark_a);
    name_to_function_object.emplace("foc_b", flow_over_cylinder_benchmark_b);
    name_to_function_object.emplace("kov_u", kov_u);
    name_to_function_object.emplace("kov_v", kov_v);
    name_to_function_object.emplace("kov_p", kov_p);
}

Orchestrator::~Orchestrator(){
}

/**
 * * Lê todas as informações do arquivo e armazena em memória.
 */
void Orchestrator::readYamlAndRecoverVariables(string yaml_filepath){
    YAML::Node config = YAML::LoadFile(yaml_filepath);
    

    // problem statement
    YAML::Node problem_definitions = config["problem"];
    this->name = problem_definitions["name"].as<string>();
    this->nu = problem_definitions["nu"].as<double>();
    this->rho = problem_definitions["rho"].as<double>();
    this->reynolds = problem_definitions["reynolds"].as<double>();

    // boundaries
    YAML::Node boundary_conditions = problem_definitions["boundaries"];
    YAML::Node u_boundaries = boundary_conditions["u"];
    for (std::size_t i = 0; i < u_boundaries.size(); ++i) {
        YAML::Node item = u_boundaries[i];

        auto t = item["type"].as<string>();
        BoundaryType aux;
        if(t == "DIRICHLET")
            aux = DIRICHLET;
        else if(t == "NEUMANN")
            aux = NEUMANN;

        bcsu.push_back(BoundaryCondition(aux, item["region"].as<string>(), name_to_function_object[item["value"].as<string>()]));
    }
    YAML::Node v_boundaries = boundary_conditions["v"];
    for (std::size_t i = 0; i < v_boundaries.size(); ++i) {
        YAML::Node item = v_boundaries[i];
        
        auto t = item["type"].as<string>();
        BoundaryType aux;
        if(t == "DIRICHLET")
            aux = DIRICHLET;
        else if(t == "NEUMANN")
            aux = NEUMANN;
        bcsv.push_back(BoundaryCondition(aux, item["region"].as<string>(), name_to_function_object[item["value"].as<string>()]));
    }
    
    YAML::Node p_boundaries = boundary_conditions["p"];
    for (std::size_t i = 0; i < p_boundaries.size(); ++i) {
        YAML::Node item = p_boundaries[i];
        
        auto t = item["type"].as<string>();
        BoundaryType aux;
        if(t == "DIRICHLET")
            aux = DIRICHLET;
        else if(t == "NEUMANN")
            aux = NEUMANN;
        bcsp.push_back(BoundaryCondition(aux, item["region"].as<string>(), name_to_function_object[item["value"].as<string>()]));
    }
    
    
    this->mshfile = config["mshfile"].as<string>();
    YAML::Node simple = config["simple"];
    YAML::Node mom = simple["momentum"];
    this->lambda_uv = mom["lambda_uv"].as<double>();
    this->non_corrections = mom["non_corrections"].as<int>();
    this->iterations_mom = mom["iterations_bicgstab"].as<int>();
    this->tolerance_mom = mom["tolerance_bicgstab"].as<double>();
    
    YAML::Node pc = simple["pressure_correction"];
    this->lambda_p = pc["lambda_p"].as<double>();
    this->iterations_pc = pc["iterations_bicgstab"].as<int>();
    this->tolerance_pc = pc["tolerance_bicgstab"].as<double>();

    this->utol = simple["utol"].as<double>();
    this->vtol = simple["vtol"].as<double>();
    this->ptol = simple["ptol"].as<double>();

    this->save_iterations = config["save_iterations"].as<bool>();
    this->exportfolder =  config["exportfolder"].as<string>();
}