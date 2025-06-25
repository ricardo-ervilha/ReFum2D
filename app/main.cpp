#include "FVMSolver.h"

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
    return 0.0;
}

pair<double,double> U(double x, double y){
    return make_pair(0.0,0.0);
}

double source(double x, double y){
    return -200*x*(1-x) - 200*y*(1-y);
}

double exact(double x, double y){
    return 100*x*(1-x)*y*(1-y);
}

int main(void){
    Mesh* m = new Mesh();
    m->read_mesh("../inputs/20x20.msh");

    BoundaryCondition* downBC = new BoundaryCondition("Dirichlet", down);
    BoundaryCondition* rightBC = new BoundaryCondition("Dirichlet", right);
    BoundaryCondition* topBC = new BoundaryCondition("Dirichlet", top);
    BoundaryCondition* leftBC = new BoundaryCondition("Dirichlet", left);
    FVMSolver* solver = new FVMSolver(m, downBC, rightBC, topBC, leftBC, gamma, rho, U, source);

    solver->assembly_A();
    // solver->print_A();

    solver->assembly_b();
    // solver->print_b();

    double tol = 1e-6;
    solver->solve_system(tol);
    solver->compute_error(exact);

    solver->save_solution("../outputs/result.vtk");

    return 0;
}