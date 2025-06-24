#ifndef PHYSICALGROUPS_H
#define PHYSICALGROUPS_H

#include "pch.h"
#include "Element.h"

class PhysicalGroup{
    private:
        vector<Element*> elements; 
        
    public:
        const int id;
        const int dimension; // 1: linha, 2: superfície
        const string name; 
    
        PhysicalGroup(int id, int dimension, string name) : id(id), dimension(dimension), name(name) {};
        PhysicalGroup() : id(), dimension(), name() {}; // construtor vazio por causa da map.
        ~PhysicalGroup() {};

        void add_element(Element* e)    {this->elements.push_back(e);};
        vector<Element*>& get_elements()                {return this->elements;};
};


#endif