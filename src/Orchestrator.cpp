#include "Orchestrator.h"
#include <yaml-cpp/yaml.h>

Orchestrator::Orchestrator(){
    name_to_function_object.emplace("one", one);
    name_to_function_object.emplace("zero", zero);
    name_to_function_object.emplace("bfs", backward_facing_step);
    name_to_function_object.emplace("foc_a", flow_over_cylinder_benchmark_a);
    name_to_function_object.emplace("foc_b", flow_over_cylinder_benchmark_b);
}

Orchestrator::~Orchestrator(){
}

/**
 * * Lê todas as informações do arquivo e armazena em memória.
 */
void Orchestrator::readYamlAndRecoverVariables(string yaml_filepath){
    YAML::Node config = YAML::LoadFile(yaml_filepath);
    
    this->problem = config["problem"].as<string>();
    this->mesh_path = config["mesh"].as<string>();
    
    YAML::Node parameters = config["parameters"];
    this->rho = parameters["density"].as<double>();
    this->nu = parameters["viscosity"].as<double>();

    YAML::Node boundary_conditions = config["boundaries"];
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
    
    auto saux = config["solver"].as<string>();
    if(saux == "transient")
        this->st = Transient;
    else if(saux == "steady")
        this->st = Steady;
    
    if(this->st == Transient){
        YAML::Node transient_parameters = config["transient_parameters"];
        this->n_steps = transient_parameters["n_steps"].as<int>();
        this->tf = transient_parameters["tf"].as<double>();
        
        YAML::Node initial_condition = transient_parameters["initial_condition"];
        
        u_ic = name_to_function_object[initial_condition["u"].as<string>()];
        v_ic = name_to_function_object[initial_condition["v"].as<string>()];
        p_ic = name_to_function_object[initial_condition["p"].as<string>()];
    }
    
    YAML::Node simple_config = config["simple"];
    this->lambda_p = simple_config["lambda_p"].as<double>();
    this->lambda_uv = simple_config["lambda_uv"].as<double>();
    this->iterations = simple_config["iterations"].as<int>();
    
    this->export_path = config["export_folder"].as<string>();
}