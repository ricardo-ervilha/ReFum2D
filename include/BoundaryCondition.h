#ifndef BOUNDARYCONDITION_H
#define BOUNDARYCONDITION_H

#include "pch.h"

enum BoundaryType {
    DIRICHLET = 1,
    NEUMANN = 2
};

class BoundaryCondition {
    private:
        string name; // O nome da boundary condition
        BoundaryType type; //Infere o tipo numérico baseado no name, que pode ser: "Dirichlet" ou "Neumann"
        double (*func)(double, double); // O valor a ser aplicado é uma função.

    public: 
        BoundaryCondition(string nameType, double (*f)(double, double)) {
            this->name = nameType;
            this->func = f;

            if(this->name == "Dirichlet")
                this->type = DIRICHLET;
            else if(this->name == "Neumann")
                this->type = NEUMANN;
            else
                throw invalid_argument("BoundaryCondition: Incorrect type of boundary condition.");
        };
        ~BoundaryCondition() {};    
        int get_type() { return type; }; 
        double apply(double x, double y) {  return func(x, y); };
};

#endif