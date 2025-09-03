#ifndef MESH_H
#define MESH_H

#include "pch.h"
class PhysicalEntity;
class Node;
class Cell;
class Edge;

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

        // variáveis para ajudar a aplicar BCs.
        double xmin, xmax;
        double ymin, ymax;

        /*Mapeia de elemento p/ o número de nós que o elemento possui*/
        map<int,int> elementTypeToNumNodes;

        /*Função para realizar o pré-processamento e obter informações úteis*/
        void pre_processing();
    public:
        Mesh();
        ~Mesh();
        void read_mesh(string filepath);

        vector<Node*>& get_nodes() { return nodes; }
        vector<Cell*>& get_cells() { return cells; }
        vector<Edge*>& get_edges() { return edges; }
        map<int, PhysicalEntity*>& getPhysicalEntities() { return physicalentities; }

        int get_nnodes() { return nodes.size(); }
        int get_ncells()  { return cells.size(); }
        int get_nedges()  { return edges.size(); }
        int get_nphysicalentities()  { return physicalentities.size(); }

        int get_xmin() {return xmin;};
        int get_xmax() {return xmax;};
        int get_ymin() {return ymin;};
        int get_ymax() {return ymax;};
};      

#endif