#include "Mesh.h"
#include "NSSolver.h"
#include "BoundaryCondition.h"
#include "utils.h"

int main(int argc, char* argv[]) {
    Mesh m;
    m.read_mesh("../inputs/q5x5.msh");

    vector<BoundaryCondition> us;
    vector<BoundaryCondition> vs;
    vector<BoundaryCondition> ps;

    us.push_back(BoundaryCondition(DIRICHLET, top_check,fone));
    us.push_back(BoundaryCondition(DIRICHLET, bottom_check, fzero));
    us.push_back(BoundaryCondition(DIRICHLET, right_check, fzero));
    us.push_back(BoundaryCondition(DIRICHLET, left_check, fzero));

    vs.push_back(BoundaryCondition(DIRICHLET, top_check,fzero));
    vs.push_back(BoundaryCondition(DIRICHLET, bottom_check, fzero));
    vs.push_back(BoundaryCondition(DIRICHLET, right_check, fzero));
    vs.push_back(BoundaryCondition(DIRICHLET, left_check, fzero));

    ps.push_back(BoundaryCondition(NEUMANN, top_check,fone));
    ps.push_back(BoundaryCondition(NEUMANN, bottom_check, fzero));
    ps.push_back(BoundaryCondition(NEUMANN, right_check, fzero));
    ps.push_back(BoundaryCondition(NEUMANN, left_check, fzero));

    NSSolver solver(&m, 1e-3, 1.0, us, vs, ps);
    
    // for(int i = 0; i < 400; i++){
    //     cout << "Calcula A_mom, b_mom_x e b_mom_y\n";
    //     solver.mom_links_and_sources(0.8);
    //     cout << "Resolve x_mom\n";
    //     solver.solve_x_mom();
    //     cout << "Resolve y_mom\n";
    //     solver.solve_y_mom();
        
    //     cout << "Calcula velocidade nas faces\n";
    //     solver.face_velocity();
    //     cout << "Calcula pressure correction\n";
    //     solver.solve_pp();
    //     cout << "Atualiza velocidades\n";
    //     solver.uv_correct();
    //     cout << "Atualiza pressão\n";
    //     solver.pres_correct(0.1);
    // }
    
    solver.export_solution("../outputs/v_vector.vtk");
}
