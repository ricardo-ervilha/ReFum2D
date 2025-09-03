#ifndef PHYSICALENTITY_H
#define PHYSICALENTITY_H

#include "pch.h"
class Cell;

class PhysicalEntity{
    private:
        vector<Cell*> cells; 
        
    public:
        const int id;
        const string name; 
    
        PhysicalEntity(int id, string name) : id(id), name(name) {};
        PhysicalEntity() : id(), name() {}; // construtor vazio por causa da map.
        ~PhysicalEntity() {};

        void add_element(Cell* e)    {this->cells.push_back(e);};
        vector<Cell*>& get_elements() {return this->cells;};
};


#endif