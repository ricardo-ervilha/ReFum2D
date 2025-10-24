#include "Orchestrator.h"
#include <yaml-cpp/yaml.h>

Orchestrator::Orchestrator(string yaml_filepath){
    this->readYamlAndRecoverVariables(yaml_filepath);
}

Orchestrator::~Orchestrator(){
    // nada por enquanto
}

