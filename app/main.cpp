#include "Mesh.h"
#include "NSSolver.h"
#include "BoundaryCondition.h"

int main(int argc, char* argv[]) {
    Mesh m;
    m.read_mesh("../inputs/q5x5.msh");
    
    auto lid = [](double x, double y) {
        return 1.0;
    };
    auto zero = [](double x, double y) {
        return 0.0;
    };

    BoundaryCondition uDown(DIRICHLET, DOWN, zero, u);
    BoundaryCondition uRight(DIRICHLET, RIGHT, zero, u);
    BoundaryCondition uTop(DIRICHLET, TOP, lid, u);
    BoundaryCondition uLeft(DIRICHLET, LEFT, zero, u);

    BoundaryCondition vDown(DIRICHLET, DOWN, zero, v);
    BoundaryCondition vRight(DIRICHLET, RIGHT, zero, v);
    BoundaryCondition vTop(DIRICHLET, TOP, zero, v);
    BoundaryCondition vLeft(DIRICHLET, LEFT, zero, v);

    BoundaryCondition pDown(NEUMANN, DOWN, zero, v);
    BoundaryCondition pRight(NEUMANN, RIGHT, zero, v);
    BoundaryCondition pTop(NEUMANN, TOP, zero, v);
    BoundaryCondition pLeft(NEUMANN, LEFT, zero, v);
    
    
    NSSolver solver(&m, 1e-2, 1.0, &uDown, &uRight, &uTop, &uLeft, &vDown, &vRight, &vTop, &vLeft, &pDown, &pRight, &pTop, &pLeft);
    
    for(int i = 0; i < 1; i++){
        cout << "Calcula A_mom, b_mom_x e b_mom_y\n";
        solver.mom_links_and_sources(0.8);
        cout << "Resolve x_mom\n";
        solver.solve_x_mom();
        cout << "Resolve y_mom\n";
        // solver.solve_y_mom();
        
        // cout << "Calcula velocidade nas faces\n";
        // solver.face_velocity();
        // cout << "Calcula pressure correction\n";
        // solver.solve_pp();
        // cout << "Atualiza velocidades\n";
        // solver.uv_correct();
        // cout << "Atualiza pressão\n";
        // solver.pres_correct(0.1);
    }
    
    solver.export_solution("../outputs/v_vector.vtk");
}
