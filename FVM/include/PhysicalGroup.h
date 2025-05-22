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
        PhysicalGroup(){
            this->dimension = 0;
            this->id = 0;
            this->name = "";
        }
        ~PhysicalGroup() {};
        void insertElementId(int id) {
            this->elementIds.push_back(id);
        }
        int getDimension()  {return this->dimension;};
        int getId()         {return this->id;};
        string getName()    {return this->name;};
        vector<int>* getElementIds()    {return &this->elementIds;};
};


#endif