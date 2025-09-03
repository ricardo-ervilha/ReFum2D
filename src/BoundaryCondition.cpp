#include "BoundaryCondition.h"
#include "Edge.h"

double BoundaryCondition::xmin;
double BoundaryCondition::xmax;
double BoundaryCondition::ymin;
double BoundaryCondition::ymax;

BoundaryLocation BoundaryCondition::find_location(Edge* e) {
    if(e->from->y == BoundaryCondition::ymin && e->to->y == BoundaryCondition::ymin){
        return DOWN;
    } else if(e->from->x == BoundaryCondition::xmax && e->to->x == BoundaryCondition::xmax){
        return RIGHT;
    } else if(e->from->y == BoundaryCondition::ymax && e->to->y == BoundaryCondition::ymax){
        return TOP;
    } else {
        return LEFT;
    }
}