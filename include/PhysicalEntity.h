#ifndef PHYSICALENTITY_H
#define PHYSICALENTITY_H

#include "pch.h"
class Cell;

class PhysicalEntity{        
    public:
        const int id;
        const string name; 
    
        PhysicalEntity(int id, string name) : id(id), name(name) {};
        PhysicalEntity() : id(), name() {}; // construtor vazio por causa da map.
        ~PhysicalEntity() {};
};


#endif