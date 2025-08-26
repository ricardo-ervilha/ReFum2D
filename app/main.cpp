#include "../include/FVMSolver.h"

double down(double x, double y) {
    return x*x + y;
}

double right(double x, double y) {
    return 2.0 * x;
}

double top(double x, double y) {
    return  x*x + y;
}

double left(double x, double y) {
    return x*x + y;
}

double gamma(double x, double y){
    return 1.0;
}

double rho(double x, double y){
    return 1.0; 
}

pair<double,double> U(double x, double y){
    return make_pair(1.0, 1.0);
}

double source(double x, double y){
    return -2.0;
}

double exact(double x, double y){
    return x*x + y; 
}

int main(void){
    Mesh* m = new Mesh();
    m->read_mesh("../inputs/100x100.msh");

    BoundaryCondition* downBC = new BoundaryCondition(DIRICHLET, DOWN, down);
    BoundaryCondition* rightBC = new BoundaryCondition(NEUMANN, RIGHT, right);
    BoundaryCondition* topBC = new BoundaryCondition(DIRICHLET, TOP, top);
    BoundaryCondition* leftBC = new BoundaryCondition(DIRICHLET, LEFT, left);
    FVMSolver* solver = new FVMSolver(m, downBC, rightBC, topBC, leftBC, gamma, rho, U, source);

    solver->diffusion(); // Calcula A e b com difusão fixa
    // solver->convection(); // Adiciona em A e b convecção

    double tol = 1e-8;
    solver->solve_system(tol); // resolve A e b, no meio do caminho computando skew e corrigindo.

    solver->compute_error(exact);
    solver->save_solution("../outputs/result.vtk");

    /* Limpando ponteiros antes declarados */
    delete m;
    delete downBC;
    delete rightBC;
    delete topBC;
    delete leftBC;
    delete solver;

    return 0;
}