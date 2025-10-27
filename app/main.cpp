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
    cout << m.get_ncells() << endl;
    
    ReFumSolver solver(&m, o.get_nu(), o.get_rho(), o.get_ubcs(), o.get_vbcs(), o.get_pbcs(), o.get_solver_type());

    // sendo transiente, é necessário passar condição inicialç
    // solver.set_initial_condition(o.get_uic(), o.get_vic(), o.get_pic());

    solver.STEADY_SIMPLE(o.get_problem(), o.get_export_path(), o.get_iterations(), o.get_lambda_uv(), o.get_lambda_p(), true);
}
