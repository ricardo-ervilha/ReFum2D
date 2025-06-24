#ifndef MESH_H
#define MESH_H

#include "pch.h"
#include "PhysicalGroup.h"
#include "Element.h"
#include "Node.h"
#include "BoundaryCondition.h"
#include "Cell.h"

// Função de hash para a parte das faces únicas.
struct CantorPairingFunction {
    // # injetora para inteiros >= 0.
    std::size_t operator()(const std::pair<int, int>& p) const noexcept {
        int a = p.first;
        int b = p.second;
        return (a + b) * (a + b + 1) / 2 + b;
    }
};

class Mesh{
    private:
        vector<Node*> nodes;
        vector<Element*> elements;
        vector<Cell*> cells;
        vector<Edge*> edges;
        map<int, PhysicalGroup*> physicalgroups;

        /*Mapeia de elemento p/ o número de nós que o elemento possui*/
        map<int,int> elementTypeToNumNodes;

        /*Função para realizar o pré-processamento e obter informações úteis*/
        void pre_processing();
    public:
        Mesh();
        ~Mesh();
        void read_mesh(string filepath);
};      

#endif