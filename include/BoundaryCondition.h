#ifndef BOUNDARYCONDITION_H
#define BOUNDARYCONDITION_H

#include <functional>
class Edge;

enum BoundaryType {
    DIRICHLET,
    NEUMANN,
    NONE
};

class BoundaryCondition {
    private:
        BoundaryType type; 
        string location;
        function<double(double, double)> func; // O valor a ser aplicado é uma função.
        static double xmin, xmax, ymin, ymax;

    public: 
        BoundaryCondition(BoundaryType type, string location, function<double(double, double)> f) {
            this->func = f;
            this->location = location;
            this->type = type;
        };
        ~BoundaryCondition() {};    
        BoundaryType get_type() { return type; }; 
        string get_location() { return location; }; 
        double apply(double x, double y) { return func(x, y); };
};

#endif