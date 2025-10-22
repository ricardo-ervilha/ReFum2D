#include "Mesh.h"
#include "NSSolver.h"
#include "BoundaryCondition.h"
#include "utils.h"

double parabolic(double x, double y){
    if(y > 0.5) // metade superior
        return -16*pow((y - 0.75), 2) + 1;
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

    m.read_mesh("../inputs/q30x30.msh");
    us.push_back(BoundaryCondition(DIRICHLET, top_check,fone));
    us.push_back(BoundaryCondition(DIRICHLET, bottom_check, fzero));
    us.push_back(BoundaryCondition(DIRICHLET, right_check, fzero));
    us.push_back(BoundaryCondition(DIRICHLET, left_check, fzero));
    
    vs.push_back(BoundaryCondition(DIRICHLET, top_check,fzero));
    vs.push_back(BoundaryCondition(DIRICHLET, bottom_check, fzero));
    vs.push_back(BoundaryCondition(DIRICHLET, right_check, fzero));
    vs.push_back(BoundaryCondition(DIRICHLET, left_check, fzero));
    
    ps.push_back(BoundaryCondition(NEUMANN, top_check,fzero));
    ps.push_back(BoundaryCondition(NEUMANN, bottom_check, fzero));
    ps.push_back(BoundaryCondition(NEUMANN, right_check, fzero));
    ps.push_back(BoundaryCondition(NEUMANN, left_check, fzero));
    NSSolver solver(&m, 1e-2, 1.0, us, vs, ps);
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
    // m.read_mesh("../inputs/backward_step_laminar_with_parabolic_function.msh");
    // us.push_back(BoundaryCondition(DIRICHLET, top_check,fzero)); // no-slip
    // us.push_back(BoundaryCondition(DIRICHLET, bottom_check, fzero)); // no-slip
    // us.push_back(BoundaryCondition(NEUMANN, right_check, fzero)); // outlet (neumann 0)
    // us.push_back(BoundaryCondition(DIRICHLET, left_check, parabolic)); // inlet u = parabolic(x,y)
    // us.push_back(BoundaryCondition(DIRICHLET, step_1_check, fzero)); // no-slip
    // us.push_back(BoundaryCondition(DIRICHLET, step_2_check, fzero)); // no-slip

    // vs.push_back(BoundaryCondition(DIRICHLET, top_check, fzero)); // no-slip
    // vs.push_back(BoundaryCondition(DIRICHLET, bottom_check, fzero)); //no-slip
    // vs.push_back(BoundaryCondition(NEUMANN, right_check, fzero)); // *outlet (neumann 0)
    // vs.push_back(BoundaryCondition(DIRICHLET, left_check, fzero)); // *inlet v = 0
    // vs.push_back(BoundaryCondition(DIRICHLET, step_1_check, fzero)); //no-slip
    // vs.push_back(BoundaryCondition(DIRICHLET, step_2_check, fzero)); //no-slip

    // ps.push_back(BoundaryCondition(NEUMANN, top_check,fzero)); //no-slip
    // ps.push_back(BoundaryCondition(NEUMANN, bottom_check, fzero)); // no-slip
    // ps.push_back(BoundaryCondition(DIRICHLET, right_check, fzero)); // *outlet (dirichlet = 0)
    // ps.push_back(BoundaryCondition(NEUMANN, left_check, fzero)); // *inlet = neumann 0
    // ps.push_back(BoundaryCondition(NEUMANN, step_1_check, fzero)); // no-slip
    // ps.push_back(BoundaryCondition(NEUMANN, step_2_check, fzero)); // no-slip
    // NSSolver solver(&m, 2e-2, 1.0, us, vs, ps);
    // & Backward facing step ----------------------------------------------------
    solver.TransientSimple();
}

// adicionar o controle do pp para true e false quando for lid e tomar cuidado com os lambdas pq o backward dependendo diverge...