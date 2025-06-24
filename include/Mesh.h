#ifndef MESH_H
#define MESH_H

#include "pch.h"
#include "PhysicalEntity.h"
#include "Node.h"
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
        vector<Cell*> cells;
        vector<Edge*> edges;
        map<int, PhysicalEntity*> physicalentities; //fluido, placa, aerofolio, etc.

        /*Mapeia de elemento p/ o número de nós que o elemento possui*/
        map<int,int> elementTypeToNumNodes;

        /*Função para realizar o pré-processamento e obter informações úteis*/
        void pre_processing();
    public:
        Mesh();
        ~Mesh();
        void read_mesh(string filepath);

        vector<Node*>& get_nodes() const { return nodes; }
        vector<Cell*>& get_cells() const { return cells; }
        vector<Edge*>& get_edges() const { return edges; }
        map<int, PhysicalEntity*>& getPhysicalEntities() const { return physicalentities; }

        int get_nnodes() const { return nodes.size(); }
        int get_ncells() const { return cells.size(); }
        int get_nedges() const { return edges.size(); }
        int get_nphysicalentities() const { return physicalentities.size(); }
};      

#endif