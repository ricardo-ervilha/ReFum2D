#ifndef BOUNDARYCONDITION_H
#define BOUNDARYCONDITION_H

#include <iostream>
#include <string>

using namespace std;

class BoundaryCondition {
    private:
        string name;
        int type;
        double value;

    public:
        BoundaryCondition(string name, double value) {
            this->name = name;
            this->value = value;

            if(this->name == "Dirichlet")
                this->type = 1;
            else if(this->name == "Neumann")
                this->type = 2;
        };
        ~BoundaryCondition() {};
        int getType() {return type;};
};

#endif