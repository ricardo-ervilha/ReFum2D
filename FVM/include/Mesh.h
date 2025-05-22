#ifndef MESH_H
#define MESH_H

#include <iostream>
#include <vector>
#include <map>
#include "Node.h"
#include "PhysicalGroup.h"
#include "Element.h"

using namespace std;

class Mesh{
    private:
        int geom_type; //Tipo da geometria (2d ou 3d)
        int ncells; //Qtd total de células
        int nfaces; //Qtd total de faces
        int nbfaces; //Qtd total de faces de contorno
        int nnodes; //Qtd total de nós/vértices
        int totalElements; //Qtd total de elementos 
        int totalPhysicalEntities; //Qtd total de grupos físicos

        map<int, Node> nodes; //int NodeID => node, contendo suas posições em x, y e z.
        map<int, PhysicalGroup> physicalGroups; //int PhysicalGroupID => grupo físico
        map<int, Element> elements; //int elementID =>  elemento da malha.

        map<int,int> elementTypeToNumNodes; // int id do elemento => int número de nós que o compõem
    public:
        Mesh();
        ~Mesh();
        void readMesh(string filepath); // Função para ler os dados da malha .msh
        void meshSummary();
};

#endif