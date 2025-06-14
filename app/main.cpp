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
    return sin(M_PI * x) * sinh(M_PI);
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
    return 0.0;
}

int main(void){
    // Declara um ponteiro para malha (Idealmente depois era passar isso pra dentro do FVMSolver... e só fornecer o path)
    Mesh* mesh = new Mesh();
    mesh->read_mesh("../inputs/mesh1.msh");

    // Opção de gerar um sumário com as estatísticas da malha
    // mesh->mesh_summary();

    BoundaryCondition down = BoundaryCondition("Dirichlet", DirichletDown);
    BoundaryCondition right = BoundaryCondition("Dirichlet", DirichletRight);
    BoundaryCondition top = BoundaryCondition("Dirichlet", DirichletTop);
    BoundaryCondition left = BoundaryCondition("Dirichlet", DirichletLeft);    
    
    /*Aplicação das condições de contorno: passar sempre em sentido anti-horário a partir do down.*/
    FVMSolver* solver = new FVMSolver(mesh, &down, &right, &top, &left, Gamma, Q);

    /* A partir daqui já temos a malha com as respectivas condições de contorno implementadas. */
}