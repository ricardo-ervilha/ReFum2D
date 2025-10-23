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

double u_velocity_cylinder(double x, double y){
    return 4 * 0.3 * y * (0.41 - y)/(0.41*0.41);
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
    
    // & Flow Over Cylinder ------------------------------------------------------
    m.read_mesh("../inputs/cylinder/flow_over_cylinder_benchmark_refined.msh");
    us.push_back(BoundaryCondition(DIRICHLET, top_check,fzero)); // no-slip
    us.push_back(BoundaryCondition(DIRICHLET, bottom_check, fzero)); // no-slip
    us.push_back(BoundaryCondition(NEUMANN, right_check, fzero)); // outlet (neumann 0)
    us.push_back(BoundaryCondition(DIRICHLET, left_check, u_velocity_cylinder)); // inlet
    us.push_back(BoundaryCondition(DIRICHLET, circle_check, fzero)); // no-slip (OBJETO)

    vs.push_back(BoundaryCondition(DIRICHLET, top_check, fzero)); // no-slip
    vs.push_back(BoundaryCondition(DIRICHLET, bottom_check, fzero)); //no-slip
    vs.push_back(BoundaryCondition(NEUMANN, right_check, fzero)); // *outlet (neumann 0)
    vs.push_back(BoundaryCondition(DIRICHLET, left_check, fzero)); // *inlet v = 0
    vs.push_back(BoundaryCondition(DIRICHLET, circle_check, fzero)); //no-slip (OBJETO)

    ps.push_back(BoundaryCondition(NEUMANN, top_check,fzero)); //no-slip
    ps.push_back(BoundaryCondition(NEUMANN, bottom_check, fzero)); // no-slip
    ps.push_back(BoundaryCondition(DIRICHLET, right_check, fzero)); // *outlet (dirichlet = 0)
    ps.push_back(BoundaryCondition(NEUMANN, left_check, fzero)); // *inlet = neumann 0
    ps.push_back(BoundaryCondition(NEUMANN, circle_check, fzero)); // no-slip (OBJETO)
    
    // Re = 20 (steady)
    double rho = 1.0, nu = 1e-3;
    NSSolver solver(&m, nu, rho, us, vs, ps);
    // & Flow Over Cylinder ------------------------------------------------------
    
    solver.TransientSimple();
}
