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
    return 1;
}

double rho(double x, double y){
    return 1.0; 
}

pair<double,double> U(double x, double y){
    return make_pair(
        0, 
        0
    );
}

double source(double x, double y){
    return 0;
}

double exact(double x, double y){
    return x*x + y; 
}

double icFunc(double x, double y){
    double r = (x-5)*(x-5) + (y-7.5)*(y-7.5);
    return exp(-0.5 * r);
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