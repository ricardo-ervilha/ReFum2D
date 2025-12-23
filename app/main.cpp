#include "Mesh.h"
#include "ReFumSolver.h"
#include "BoundaryCondition.h"
#include "utils.h"
#include "Orchestrator.h"
#include "PhysicalEntity.h"
#include "Cell.h"

int main(int argc, char* argv[]) {
    Orchestrator o;
    o.readYamlAndRecoverVariables(argv[1]);

    // instancia malha e realiza pré-processamento
    Mesh m;
    m.read_mesh(o.get_mshfile());
    
    // informações
    cout << "=======*===========*===========*===========*===========*===========*====\n";
    cout << "# Nome do Problema: " << o.get_name() << endl;
    cout << "# Densidade: " << o.get_rho() << endl;
    cout << "# Viscosidade: " << o.get_nu() << endl;
    cout << "# Reynolds: " << o.get_reynolds() << endl;
    cout << "# Valor do lambda_uv: " << o.get_lambda_uv() << endl;
    cout << "# Valor do lambda_p: " << o.get_lambda_p() << endl;
    cout << "# Será salvado as iterações ? " << boolalpha << o.get_save_iterations() << endl;
    cout << "> Número de correções de não ortogonalidade: " << o.get_non_corrections() << endl;
    cout << "> Número de iterações BiCGStab Momentum: " << o.get_iterations_mom() << endl;
    cout << "> Número de iterações BiCGStab Pressure Correction: " << o.get_iterations_pc() << endl;
    cout << "> Tolerância BiCGStab Momentum: " << o.get_tolerance_mom() << endl;
    cout << "> Tolerância BiCGStab Pressure Correction: " << o.get_tolerance_pc() << endl;
    cout << "=======*===========*===========*===========*===========*===========*====\n";
    
    ReFumSolver solver(&m, o.get_nu(), o.get_rho(), o.get_bcsu(), o.get_bcsv(), o.get_bcsp(), Steady);
    
    solver.SIMPLE(o.get_name(), o.get_reynolds(), o.get_exportfolder(), o.get_save_iterations(), o.get_lambda_uv(), o.get_lambda_p(), o.get_non_corrections(), o.get_iterations_mom(), o.get_tolerance_mom(), o.get_iterations_pc(), o.get_tolerance_pc(), o.get_utol(), o.get_vtol(), o.get_ptol());
}
