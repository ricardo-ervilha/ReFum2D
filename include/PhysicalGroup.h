#ifndef PHYSICALGROUPS_H
#define PHYSICALGROUPS_H

#include <iostream>

using namespace std;

class PhysicalGroup{
    private:
        int dimension; // 1: linha, 2: superfície, 3: volume
        int id; // identificador numérico
        string name; //nome para referir a condição de contorno
        vector<int> elementIds;
    
    public:
        PhysicalGroup(int dimension, int id, string name) {
            this->dimension = dimension;
            this->id = id;
            this->name = name;
        };
        PhysicalGroup(){ //construtor default para o map
            this->dimension = 0;
            this->id = 0;
            this->name = "";
        }
        ~PhysicalGroup() {};
        
        int get_dimension()  {return this->dimension;};
        int get_id()         {return this->id;};
        string get_name()    {return this->name;};
        vector<int>& get_element_ids()    {return this->elementIds;};
        
        void insert_element_id(int id) {this->elementIds.push_back(id);}
};


#endif