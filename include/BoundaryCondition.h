#ifndef BOUNDARYCONDITION_H
#define BOUNDARYCONDITION_H

class Edge;

enum BoundaryType {
    DIRICHLET,
    NEUMANN
};

enum BoundaryLocation {
    DOWN,
    RIGHT, 
    TOP,
    LEFT
};

class BoundaryCondition {
    private:
        BoundaryType type; 
        BoundaryLocation location;
        double (*func)(double, double); // O valor a ser aplicado é uma função.
        static double xmin, xmax, ymin, ymax;

    public: 
        BoundaryCondition(BoundaryType type, BoundaryLocation location, double (*f)(double, double)) {
            this->func = f;
            this->location = location;
            this->type = type;
        };
        ~BoundaryCondition() {};    
        BoundaryType get_type() { return type; }; 
        BoundaryLocation get_location() { return location; }; 
        double apply(double x, double y) {  return func(x, y); };

        static void apply_mins_maxs(double xmin, double xmax, double ymin, double ymax){
            BoundaryCondition::xmin = xmin;
            BoundaryCondition::xmax = xmax;
            BoundaryCondition::ymin = ymin;
            BoundaryCondition::ymax = ymax;
        }

        static BoundaryLocation find_location(Edge* e);
};

#endif