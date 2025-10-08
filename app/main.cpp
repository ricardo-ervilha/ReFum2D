#include "Mesh.h"
#include "NSSolver.h"
#include "BoundaryCondition.h"
#include <cmath>

// double u_func(double x, double y) {
//     double delta = 0.5 - sqrt(0.25 + 4 * M_PI * M_PI);
//     return 1.0 - std::exp(delta * x) * std::cos(2 * M_PI * y);
// }

// double v_func(double x, double y) {
//     double delta = 0.5 - sqrt(0.25 + 4 * M_PI * M_PI);
//     return (delta / (2 * M_PI)) * std::exp(delta * x) * std::sin(2 * M_PI * y);
// }

// double p_func(double x, double y) {
//     double delta = 0.5 - sqrt(0.25 + 4 * M_PI * M_PI);
//     return 0.5 * std::exp(2 * delta * x);
// }

double lid(double x, double y){
    return 1.0;
}

double zero(double x, double y){
    return 0.0;
}

int main(int argc, char* argv[]) {
    Mesh m;
    m.read_mesh("../inputs/q30x30.msh");

    // BoundaryCondition uDown(DIRICHLET, DOWN, u_func, u);
    // BoundaryCondition uRight(DIRICHLET, RIGHT, u_func, u);
    // BoundaryCondition uTop(DIRICHLET, TOP, u_func, u);
    // BoundaryCondition uLeft(DIRICHLET, LEFT, u_func, u);

    // BoundaryCondition vDown(DIRICHLET, DOWN, v_func, v);
    // BoundaryCondition vRight(DIRICHLET, RIGHT, v_func, v);
    // BoundaryCondition vTop(DIRICHLET, TOP, v_func, v);
    // BoundaryCondition vLeft(DIRICHLET, LEFT, v_func, v);

    // BoundaryCondition pDown(DIRICHLET, DOWN, p_func, p);
    // BoundaryCondition pRight(DIRICHLET, RIGHT, p_func, p);
    // BoundaryCondition pTop(DIRICHLET, TOP, p_func, p);
    // BoundaryCondition pLeft(DIRICHLET, LEFT, p_func, p);
    
    BoundaryCondition uDown(DIRICHLET, DOWN, zero, u);
    BoundaryCondition uRight(DIRICHLET, RIGHT, zero, u);
    BoundaryCondition uTop(DIRICHLET, TOP, lid, u);
    BoundaryCondition uLeft(DIRICHLET, LEFT, zero, u);

    BoundaryCondition vDown(DIRICHLET, DOWN, zero, v);
    BoundaryCondition vRight(DIRICHLET, RIGHT, zero, v);
    BoundaryCondition vTop(DIRICHLET, TOP, zero, v);
    BoundaryCondition vLeft(DIRICHLET, LEFT, zero, v);

    BoundaryCondition pDown(NEUMANN, DOWN, zero, p);
    BoundaryCondition pRight(NEUMANN, RIGHT, zero, p);
    BoundaryCondition pTop(NEUMANN, TOP, zero, p);
    BoundaryCondition pLeft(NEUMANN, LEFT, zero, p);
    
    
    NSSolver solver(&m, 1e-2, 1.0, &uDown, &uRight, &uTop, &uLeft, &vDown, &vRight, &vTop, &vLeft, &pDown, &pRight, &pTop, &pLeft);
    
    for(int i = 0; i < 100; i++){
        cout << "Calcula A_mom, b_mom_x e b_mom_y\n";
        solver.mom_links_and_sources(0.8);
        cout << "Resolve x_mom\n";
        solver.solve_x_mom();
        cout << "Resolve y_mom\n";
        // solver.solve_y_mom();
        
        // cout << "Calcula velocidade nas faces\n";
        solver.face_velocity();
        cout << "Calcula pressure correction\n";
        solver.solve_pp();
        cout << "Atualiza velocidades\n";
        solver.uv_correct();
        cout << "Atualiza pressão\n";
        solver.pres_correct(0.1);
    }
    
    solver.export_solution("../outputs/v_vector.vtk");
}
