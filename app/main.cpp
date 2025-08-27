#include "../include/FVMSolver.h"

double down(double x, double y) {
    return 0.0;
}

double right(double x, double y) {
    return 1.0;
}

double top(double x, double y) {
    return  0.0;
}

double left(double x, double y) {
    return 0.0;
}

double gamma(double x, double y){
    return 0.002;
}

double rho(double x, double y){
    return 1; 
}

pair<double,double> U(double x, double y){
    return make_pair(
        - sin(M_PI * x) * cos(M_PI * y),
        cos(M_PI * x) * sin(M_PI * y)
    );
}

double source(double x, double y){
    return 0;
}

double exact(double x, double y){
    return x*x + y;  // ignorar
}

double icFunc(double x, double y){
    return 0;
}

int main(void){
    Mesh* m = new Mesh();
    m->read_mesh("../inputs/quadStretchRefined.msh");

    BoundaryCondition* downBC = new BoundaryCondition(NEUMANN, DOWN, down);
    BoundaryCondition* rightBC = new BoundaryCondition(DIRICHLET, RIGHT, right);
    BoundaryCondition* topBC = new BoundaryCondition(NEUMANN, TOP, top);
    BoundaryCondition* leftBC = new BoundaryCondition(DIRICHLET, LEFT, left);
    FVMSolver* solver = new FVMSolver(m, downBC, rightBC, topBC, leftBC, gamma, rho, U, source);

    solver->set_initial_condition(icFunc);
    solver->TransientSolver(0, 20, 50, "../outputs/gaussian_pulse");
    
    solver->diffusion();
    solver->convection();
    solver->SteadySolver(1e-8);
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