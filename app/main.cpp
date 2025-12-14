#include "Mesh.h"
#include "ReFumSolver.h"
#include "BoundaryCondition.h"
#include "utils.h"
#include "Orchestrator.h"
#include "PhysicalEntity.h"
#include "Cell.h"

int main(int argc, char* argv[]) {
    Orchestrator o;
    o.readYamlAndRecoverVariables("../app/backward_facing_step.yaml");

    Mesh m;
    m.read_mesh(o.get_mesh_path());
    
    cout << "# Nome do Problema: " << o.get_problem() << endl;
    if(o.get_solver_type() == Transient)
        cout << "# Tipo do problema: Transiente\n";
    else
        cout << "# Tipo do problema: Estacionário\n";
    cout << "# Densidade: " << o.get_rho() << endl;
    cout << "# Viscosidade: " << o.get_nu() << endl;

    ReFumSolver solver(&m, o.get_nu(), o.get_rho(), o.get_ubcs(), o.get_vbcs(), o.get_pbcs(), o.get_solver_type());

    if(o.get_solver_type() == Transient){
        solver.set_initial_condition(o.get_uic(), o.get_vic(), o.get_pic());
        solver.TRANSIENTE_SIMPLE(o.get_problem(), o.get_export_path(), o.get_iterations(), o.get_lambda_uv(), o.get_lambda_p(), o.get_n_steps(), o.get_tf(), false);
    }
    else{
        solver.STEADY_SIMPLE(o.get_problem(), o.get_export_path(), o.get_iterations(), o.get_lambda_uv(), o.get_lambda_p(), false);
    }    
    
}
