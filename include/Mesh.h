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

        vector<tuple<int, int, int>> boundarySegments; 

        // variáveis para ajudar a aplicar BCs.
        double xmin, xmax;
        double ymin, ymax;

        /*Mapeia de elemento p/ o número de nós que o elemento possui*/
        map<int,int> elementTypeToNumNodes;

        /*Função para realizar o pré-processamento e obter informações úteis*/
        void pre_processing();

        // & Construção de ARRAY's para melhorar o acesso aos valores (aproveitar localidade.)
        vector<vector<int>> cell_nsign;
        vector<double> cell_area;
        vector<pair<double,double>> cell_centroid;
        vector<vector<int>> cell_faces_ids;
        
        vector<double> face_length;
        vector<double> face_df;
        vector<pair<double,double>> face_middle;
        vector<pair<double,double>> face_normal;
        vector<pair<int,int>> face_lftc;
        vector<bool> face_boundary_face;

    public:
        Mesh();
        ~Mesh();
        void read_mesh(string filepath);
        
        const vector<vector<int>>& getCellNsign() const { return cell_nsign; }
        const vector<double>& getCellArea() const { return cell_area; }
        const vector<pair<double,double>>& getCellCentroid() const { return cell_centroid; }
        const vector<double>& getFaceLength() const { return face_length; }
        const vector<double>& getFaceDf() const { return face_df; }
        const vector<pair<double,double>>& getFaceMiddle() const { return face_middle; }
        const vector<pair<double,double>>& getFaceNormal() const { return face_normal; }
        const vector<pair<int,int>>& getFaceLftc() const { return face_lftc; }
        const vector<bool> getIsBoundaryFace() const {return face_boundary_face;};
        const vector<vector<int>> getFacesFromCell() const {return cell_faces_ids;};


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