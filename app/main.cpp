#include "../include/FVMSolver.h"

double down(double x, double y) {
    return 0.0;
}

double right(double x, double y) {
    return 0.0;
}

double top(double x, double y) {
    return  0.0;
}

double left(double x, double y) {
    return 0.0;
}

double gamma(double x, double y){
    return 1e-8;
}

double rho(double x, double y){
    return 1; 
}

pair<double,double> U(double x, double y){
    return make_pair(
       -y + 5.0,
       x - 5.0
    );
}

double source(double x, double y){
    return 0; 
}

double exact(double x, double y){
    return exp(x+y);  // ignorar
}

double icFunc(double x, double y){
    double r = (x-5)*(x-5) + (y-7.5)*(y-7.5);
    return exp(-0.5*r);
}

int main(void){
    Mesh* m = new Mesh();
    m->read_mesh("../inputs/refinedPulseGaussian.msh");

    BoundaryCondition* downBC = new BoundaryCondition(DIRICHLET, DOWN, down);
    BoundaryCondition* rightBC = new BoundaryCondition(DIRICHLET, RIGHT, right);
    BoundaryCondition* topBC = new BoundaryCondition(DIRICHLET, TOP, top);
    BoundaryCondition* leftBC = new BoundaryCondition(DIRICHLET, LEFT, left);
    FVMSolver* solver = new FVMSolver(m, downBC, rightBC, topBC, leftBC, gamma, rho, U, source);

    solver->set_initial_condition(icFunc);
    solver->TransientSolver(0, 2*M_PI, 500, "../outputs/gaussian_pulse");

    // solver->diffusion();
    // solver->convection();
    // solver->SteadySolver(1e-8);
    // solver->save_solution("../outputs/result.vtk");
    // solver->compute_error(exact);

    /* Limpando ponteiros antes declarados */
    delete m;
    delete downBC;
    delete rightBC;
    delete topBC;
    delete leftBC;
    delete solver;

    return 0;
}