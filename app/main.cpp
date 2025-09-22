#include "ConvectionDiffusionBenchmark.h"
#include "FVMSolver.h"
#include "Mesh.h"
#include "Diffusion.h"
#include "Convection.h"
#include "Source.h"
#include "GradientReconstruction.h"

// DIFUSÃO-CONVECÇÃO
int main(int argc, char* argv[]) {

    // argv[1]: caminho até o .msh
    // argv[2]: caminho para o .vtk
    // argv[3]: sufixo para a metrica q será salva na pasta python

    Mesh m;
    m.read_mesh(argv[1]);

    BoundaryCondition downBC(DIRICHLET, DOWN, down);
    BoundaryCondition rightBC(DIRICHLET, RIGHT, right);
    BoundaryCondition topBC(DIRICHLET, TOP, top);
    BoundaryCondition leftBC(DIRICHLET, LEFT, left);

    FVMSolver solver(&m, &downBC, &rightBC, &topBC, &leftBC);

    Diffusion d(&solver, gamma);
    Source s(&solver, source);
    GradientReconstruction g(&solver);
    Convection c(&solver, rho, U);

    d.assembleCoefficients(); // preenche os coeficientes fixos.
    s.assemblyCoefficients(); // prenche os coeficientes fixos.
    c.assembleCoefficients(); // preenche com os coeficientes fixos

    double tol = 1e-8;
    double diff = 1;
    int iter = 0;

    arma::vec aux = solver.b_aux; // pega valor antigo, tudo zero
    while(diff > tol){
        // reseta a contribuição explícita em b_aux
        solver.resetExplicit();
        
        // obtém os gradientes
        g.reconstruct_gradients();

        // a partir dos gradientes, computa difusão cruzada
        d.cross_diffusion();

        // a partir dos gradientes, computa linear upwind
        c.linear_upwind();

        // adiciona a contribuição de b em b_aux
        solver.addExplicitContribution();

        // resolve o sistema
        solver.solveSystem();

        diff = arma::norm(solver.b_aux - aux, "inf"); // computa norma infinito
        iter += 1;
        aux = solver.b_aux; // guarda o b_aux da iteração anterior.
    }

    cout << "Número de iterações para convergência: " << iter << endl;
    solver.compute_error(exact, argv[3]);
    solver.export_solution(argv[2]);

    return 0;
}

// DIFUSÃO
// int main(int argc, char* argv[]) {

//     // argv[1]: caminho até o .msh
//     // argv[2]: caminho para o .vtk
//     // argv[3]: sufixo para a metrica q será salva na pasta python

//     Mesh m;
//     m.read_mesh(argv[1]);

//     BoundaryCondition downBC(DIRICHLET, DOWN, down);
//     BoundaryCondition rightBC(DIRICHLET, RIGHT, right);
//     BoundaryCondition topBC(DIRICHLET, TOP, top);
//     BoundaryCondition leftBC(DIRICHLET, LEFT, left);

//     FVMSolver solver(&m, &downBC, &rightBC, &topBC, &leftBC);

//     Diffusion d(&solver, gamma);
//     Source s(&solver, source);
//     GradientReconstruction g(&solver);

//     d.assembleCoefficients(); // preenche os coeficientes fixos.
//     s.assemblyCoefficients(); // prenche os coeficientes fixos.

//     double tol = 1e-8;
//     double diff = 1;
//     int iter = 0;

//     arma::vec aux = solver.b_aux; // pega valor antigo, tudo zero
//     while(diff > tol){
//         // reseta a contribuição explícita em b_aux
//         solver.resetExplicit();
        
//         // obtém os gradientes
//         g.reconstruct_gradients();

//         // a partir dos gradientes, computa difusão cruzada
//         d.cross_diffusion();

//         // adiciona a contribuição de b em b_aux
//         solver.addExplicitContribution();

//         // resolve o sistema
//         solver.solveSystem();

//         diff = arma::norm(solver.b_aux - aux, "inf"); // computa norma infinito
//         iter += 1;
//         aux = solver.b_aux; // guarda o b_aux da iteração anterior.
//     }

//     cout << "Número de iterações para convergência: " << iter << endl;
//     solver.compute_error(exact, argv[3]);
//     solver.export_solution(argv[2]);

//     return 0;
// }
