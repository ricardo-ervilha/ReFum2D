#ifndef BOUNDARYCONDITION_H
#define BOUNDARYCONDITION_H

class Edge;

enum BoundaryType {
    DIRICHLET,
    NEUMANN,
    NONE
};

class BoundaryCondition {
    private:
        BoundaryType type; 
        bool (*location)(double, double);
        double (*func)(double, double); // O valor a ser aplicado é uma função.
        static double xmin, xmax, ymin, ymax;

    public: 
        BoundaryCondition(BoundaryType type, bool (*location)(double, double), double (*f)(double, double)) {
            this->func = f;
            this->location = location;
            this->type = type;
        };
        ~BoundaryCondition() {};    
        BoundaryType get_type() { return type; }; 
        bool get_location(double x, double y) { return location(x,y); }; 
        double apply(double x, double y) {  return func(x, y); };
};

#endif