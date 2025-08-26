#include "../include/FVMSolver.h"

double down(double x, double y) {
    return 0.0;
}

double right(double x, double y) {
    return 0.0;
}

double top(double x, double y) {
    return 0.0;
}

double left(double x, double y) {
    return 0.0;
}

double gamma(double x, double y){
    return 1.0;
}

double rho(double x, double y){
    return 1.0; // NN usado
}

pair<double,double> U(double x, double y){
    return make_pair(1.0, 1.0); // NN usado
}

double source(double x, double y){
    return -200*x*(1-x) - 200*y*(1-y);
}

double exact(double x, double y){
    return 100*x*(1-x)*y*(1-y); 
}

int main(void){
    Mesh* m = new Mesh();
    m->read_mesh("../inputs/540tri.msh");

    BoundaryCondition* downBC = new BoundaryCondition(DIRICHLET, DOWN, down);
    BoundaryCondition* rightBC = new BoundaryCondition(DIRICHLET, RIGHT, right);
    BoundaryCondition* topBC = new BoundaryCondition(DIRICHLET, TOP, top);
    BoundaryCondition* leftBC = new BoundaryCondition(DIRICHLET, LEFT, left);
    FVMSolver* solver = new FVMSolver(m, downBC, rightBC, topBC, leftBC, gamma, rho, U, source);

    solver->diffusion(); // Calcula A e b com difusão fixa
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