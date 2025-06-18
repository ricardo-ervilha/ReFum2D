#include "../include/Mesh.h"
#include "../include/FVMSolver.h"
#include "../include/BoundaryCondition.h"

/*==================================================================================================================*/
/*Condições de Contorno*/
double DirichletDown(double x, double y){
    return 0;
}

double DirichletRight(double x, double y){
    return 0;
}

double DirichletTop(double x, double y){
    return 0; // PROBLEMA 1
    // return sin(M_PI * x) * sinh(M_PI); //PROBLEMA 2
}

double DirichletLeft(double x, double y){
    return 0;
}
/*==================================================================================================================*/
/*Constante de difusão*/
double Gamma(double x, double y){
    return 1.0;
}
/*==================================================================================================================*/
double Q(double x, double y){
    return -200*x*(1 - x) - 200*y*(1 - y); //PROBLEMA 1
    // return 0.0; //PROBLEMA 2
}
/*==================================================================================================================*/
double exact_solution(double x, double y){
    return 100*x*(1-x)*y*(1-y); //PROBLEMA 1
}

int main(void){
    BoundaryCondition down = BoundaryCondition("Dirichlet", DirichletDown);
    BoundaryCondition right = BoundaryCondition("Dirichlet", DirichletRight);
    BoundaryCondition top = BoundaryCondition("Dirichlet", DirichletTop);
    BoundaryCondition left = BoundaryCondition("Dirichlet", DirichletLeft);    
    
    // /*Aplicação das condições de contorno: passar sempre em sentido anti-horário a partir do down.*/
    FVMSolver* solver = new FVMSolver("../inputs/6x6.msh", &down, &right, &top, &left, Gamma, Q);
    solver->summary();
    solver->assembly_A();
    solver->assembly_b();
    solver->print_A();
    solver->print_b();

    double tolerance = 1e-4;
    solver->iterative_solver(tolerance);
    solver->save_solution("../outputs/result.vtk");
    cout << endl;
    cout << "ERRO: " << solver->compute_error(exact_solution) << endl;
}