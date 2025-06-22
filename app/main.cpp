#include "FVMSolver.h"

/*==================================================================================================================*/
/*Condições de Contorno*/
double DirichletDown(double x, double y){
    return 0;
}

double DirichletRight(double x, double y){
    return 0;
}

double DirichletTop(double x, double y){
    return 0;
}

double DirichletLeft(double x, double y){
    return 0;
}
/*==================================================================================================================*/
/*Constante de difusão*/
double Gamma(double x, double y){
    return 1;
}
/*==================================================================================================================*/
/*densidade do meio*/
double Rho(double x, double y){
    return 0;
}
/*===========================================================================================================*/
//vetor velocidade do meio
pair<double,double> U(double x, double y){
    return make_pair(0, 0);
}
/*===========================================================================================================*/
/*Termo fonte*/
double Q(double x, double y){
    return -200*x*(1-x) - 200*y*(1-y);
}

int main(void){
    BoundaryCondition down = BoundaryCondition("Dirichlet", DirichletDown);
    BoundaryCondition right = BoundaryCondition("Dirichlet", DirichletRight);
    BoundaryCondition top = BoundaryCondition("Dirichlet", DirichletTop);
    BoundaryCondition left = BoundaryCondition("Dirichlet", DirichletLeft);    
    
    // /*Aplicação das condições de contorno: passar sempre em sentido anti-horário a partir do down.*/
    FVMSolver* solver = new FVMSolver("../inputs/q30x30.msh", &down, &right, &top, &left, Gamma, Rho, U, Q);

    /*Solução...*/
    solver->assembly_A();
    solver->assembly_b();

    solver->solve_system();
    solver->save_solution("../outputs/result.vtk");

    // solver->compute_error(exact_solution);
}