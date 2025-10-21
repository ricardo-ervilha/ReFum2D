#include "Mesh.h"
#include "NSSolver.h"
#include "BoundaryCondition.h"
#include "utils.h"

double parabolic(double x, double y){
    if(y > 1.0) // metade superior
        return 1.0;
    else
        return 0.0;
}

int main(int argc, char* argv[]) {
    Mesh m;
    
    vector<BoundaryCondition> us;
    vector<BoundaryCondition> vs;
    vector<BoundaryCondition> ps;
    
    // & Lid Driven Cavity Flow --------------------------------------------------
    /*
                    (0,1)           u_lid = 1.0, v = 0.0 | ∇p = 0               (1,1)
                    *------------------------------------------------------------*
                    |                                                            |
                    |                                                            |
                    |                                                            |
                    |                                                            |
                    |                                                            |
                    |                                                            |
u = v = 0 | ∇p = 0  |                                                            | u = v = 0 | ∇p = 0
                    |                                                            |
                    |                                                            |
                    |                                                            |
                    |                                                            |
                    |                                                            |
                    |                                                            |
                    |                                                            |
                    |                                                            |
                    |                                                            |
                    |                                                            |
              (0,0) *------------------------------------------------------------* (1,0)
                                        u = v = 0 | ∇p = 0
    */

    // m.read_mesh("../inputs/q30x30.msh");
    // us.push_back(BoundaryCondition(DIRICHLET, top_check,fone));
    // us.push_back(BoundaryCondition(DIRICHLET, bottom_check, fzero));
    // us.push_back(BoundaryCondition(DIRICHLET, right_check, fzero));
    // us.push_back(BoundaryCondition(DIRICHLET, left_check, fzero));
    
    // vs.push_back(BoundaryCondition(DIRICHLET, top_check,fzero));
    // vs.push_back(BoundaryCondition(DIRICHLET, bottom_check, fzero));
    // vs.push_back(BoundaryCondition(DIRICHLET, right_check, fzero));
    // vs.push_back(BoundaryCondition(DIRICHLET, left_check, fzero));
    
    // ps.push_back(BoundaryCondition(NEUMANN, top_check,fzero));
    // ps.push_back(BoundaryCondition(NEUMANN, bottom_check, fzero));
    // ps.push_back(BoundaryCondition(NEUMANN, right_check, fzero));
    // ps.push_back(BoundaryCondition(NEUMANN, left_check, fzero));
    // NSSolver solver(&m, 1e-2, 1.0, us, vs, ps);
    // & Lid Driven Cavity Flow --------------------------------------------------
    
    // & Backward facing step ----------------------------------------------------
    /*
                    (0,2)                        wall: u = 0.0, v = 0.0 | ∇p = 0                            (10,2)
                    *-----------------------------------------------------------------------------------------*
                    |                                                                                        |
                    |                                                                                        |
                    |                                                                                        |
    inflow(y > 1)   |                                                                                        | outflow
u=1, v = 0 | ∇p = 0 |                                                                                        | ∇U = 0 | ∇p = 0
                    |                                                                                        |
                    |                                                                                        |
                    |                                                                                        |
                    +-------------------+                                                                    |
                    |                   |                                                                    |
                    |                   |                                                                    |
                    |                   |                                                                    |
              (0,0) *-----------------------------------------------------------------------------------------* (10,0)
                                                      wall: u = v = 0 | ∇p = 0
    */
    m.read_mesh("../inputs/backward_complete_40_80.msh");
    us.push_back(BoundaryCondition(DIRICHLET, top_check,fzero)); // no-slip
    us.push_back(BoundaryCondition(DIRICHLET, bottom_check, fzero)); // no-slip
    us.push_back(BoundaryCondition(NEUMANN, right_check, fzero)); // outlet (neumann 0)
    us.push_back(BoundaryCondition(DIRICHLET, left_check, parabolic)); // inlet u = parabolic(x,y)

    vs.push_back(BoundaryCondition(DIRICHLET, top_check, fzero)); // no-slip
    vs.push_back(BoundaryCondition(DIRICHLET, bottom_check, fzero)); //no-slip
    vs.push_back(BoundaryCondition(NEUMANN, right_check, fzero)); // *outlet (neumann 0)
    vs.push_back(BoundaryCondition(DIRICHLET, left_check, fzero)); // *inlet v = 0

    ps.push_back(BoundaryCondition(NEUMANN, top_check,fzero)); //no-slip
    ps.push_back(BoundaryCondition(NEUMANN, bottom_check, fzero)); // no-slip
    ps.push_back(BoundaryCondition(DIRICHLET, right_check, fzero)); // *outlet (dirichlet = 0)
    ps.push_back(BoundaryCondition(NEUMANN, left_check, fzero)); // *inlet = neumann 0
    NSSolver solver(&m, 2e-2, 1.0, us, vs, ps);
    // & Backward facing step ----------------------------------------------------

    
    for(int i = 0; i < 20; i++){
        cout << "# Calculando A_mom, b_mom_x e b_mom_y\n";
        solver.mom_links_and_sources(0.6);
        cout << "Resolvendo para encontrar uc\n";
        solver.solve_x_mom();
        cout << "Resolvendo para encontrar uv\n";
        solver.solve_y_mom();
        
        cout << "# Calculando velocidade nas faces {u_f e v_f}\n";
        solver.face_velocity();
        cout << "# Calculando correção na pressão (p')\n";
        solver.solve_pp();
        cout << "# Atualiza velocidades...\n";
        solver.uv_correct();
        cout << "# Atualiza pressão....\n";
        solver.pres_correct(0.3);
    }
    
    solver.export_solution("../outputs/v_vector.vtk");
}
