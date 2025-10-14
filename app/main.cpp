#include "Mesh.h"
#include "NSSolver.h"

int main(int argc, char* argv[]) {
    Mesh m;
    m.read_mesh("../inputs/regular/quad_100x100.msh");

    NSSolver solver(&m, 1e-3, 1.0);
    
    for(int i = 0; i < 1; i++){
            cout << "Calcula A_mom, b_mom_x e b_mom_y\n";
            solver.mom_links_and_sources(0.8);
            cout << "Resolve x_mom\n";
            // solver.solve_x_mom();
            // cout << "Resolve y_mom\n";
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
    
    // solver.export_solution("../outputs/v_vector.vtk");
}
