#ifndef MESH_H
#define MESH_H

#include <iostream>
#include <vector>
#include <map>
#include <unordered_set>
#include <unordered_map>
#include <tuple>
#include <math.h>
#include "Node.h"
#include "PhysicalGroup.h"
#include "Element.h"
#include <utility>

using namespace std;

/*Função hash para uso no unordered set posteriormente...*/
struct CantorPairingFunction {
    // # injetora para inteiros >= 0.
    std::size_t operator()(const std::pair<int, int>& p) const noexcept {
        int a = p.first;
        int b = p.second;
        return (a + b) * (a + b + 1) / 2 + b;
    }
};

// Estrutura para conseguir manter a ordem original depois...
struct Edge {
    int from, to; // mantém a ordem original

    Edge(int a, int b) : from(a), to(b) {}

    pair<int, int> asKey() const {
        return (from < to) ? make_pair(from, to) : make_pair(to, from);
    }
};

class Mesh{
    private:
        /*EDs para os vértices/nós*/
        vector<Node> nodes; // guarda o tad Nó com suas informações...
        int nnodes; // registra quantidade total de nós

        /*Elements*/
        vector<Element> elements; // guarda o tad Elemeno (células e outros que vem junto do gmsh)
        int totalElements; // registra quantidade total de elementos

        /*Grupos físicos*/
        map<int, PhysicalGroup> physicalGroups; // mapeamento: int PhysicalGroupID => Grupo Físico TAD
        int totalPhysicalGroups; // registra quantidade total de grupos físicos

        /*Células*/
        int ncells; //registra quantidade total de células
        vector<int> cells; //guarda o id das células.
        vector<double> volumes; //volumes de cada célula

        /*Faces*/
        int nfaces; // registra quantidade total de faces
        vector<pair<int, int>> faces; //Vetor de faces. Ordem das faces está imbutida da indexação do vetor
        vector<pair<int,int>> link_face_to_cell; // para a i-ésima face, o pair contém as células que compartilham a face.

        /*Geral*/
        int geom_type; // Tipo da geometria, se é 2d ou 3d
        map<int,int> elementTypeToNumNodes; // int id do elemento => int número de nós que o compõem
        int offset; // utilizado para mapear de célula global para local
        
        /*Faces de contorno*/
        int nbfaces; // registra quantidad etotal de faces de contorno
        vector<int> link_face_to_bface; //vetor que acompanha a indexação das faces, e em cada entrada correspondendo a id face global, diz qual id da bface, ou retorna -1 caso não seja boundary face.
        vector<int> link_bface_to_face; // vetor do tamanho de faces de contorno. Para cada indice do vetor, diz qual o indice global da boundary face.
        unordered_map<pair<int,int>, int, CantorPairingFunction> facesPairToKey; // dado o par de chaves retorna o índice => útil na hora das faces de contorno

        /*Centroides*/
        vector<Node> centroids; // acompanha a indexação das células

        /*Normais*/
        vector<tuple<double, double>> normals; // normais de cada face. Acompanha a indexação de faces
        
        /*Área da Face (3d) ou Comprimento da face (2d)*/
        vector<double> faceAreas; // acompanha a indexação do faces.
        
        /*Ponto médio das faces*/
        vector<Node> faceMiddlePoints; //armazena o ponto central das faces. 

        /*Funções privadas*/
        void pre_processing(); //Calcula valores de diversas coisas que não foram feitos no readMesh

        /*Calcula da projeção do vetor que une dois centroides na direção normal.*/
        vector<double> deltafs; // acompanha a indexação da face 

    public:
        Mesh();
        ~Mesh();

        // Função para ler os dados da malha .msh, passar o nome do arquivo e caminho até ele.
        void read_mesh(string filepath); 
        
        void mesh_summary();

        void print_info_geral();
        void print_info_nodes();
        void print_info_elements();
        void print_info_physical_groups();
        void print_info_cells();
        void print_info_faces();
        void print_info_centroids();
        void print_info_link_face_to_cell();
        void print_info_normal_signs();
        void print_info_normal_values();
        void print_info_volumes();
        void print_info_link_face_to_bface();
        void print_info_link_bface_to_face();
        void print_deltafs();
        void print_info_distance_node_to_centroids();

        int get_num_cells()      {return this->ncells;};
        int get_num_elements()   {return this->totalElements;};
        int get_global_cell_id(int localId) {return this->cells[localId];};
        int get_link_face_to_bface(int id) {return this->link_face_to_bface[id];};
        int get_link_bface_to_face(int id) {return this->link_bface_to_face[id];};
        double get_face_area(int id) {return this->faceAreas[id];};
        double get_deltaf(int id) {return this->deltafs[id];};
        int get_num_boundary_faces() {return this->nbfaces;};
        double get_volume(int id) {return this->volumes[id];};
        int get_num_nodes() {return this->nnodes;};
        int get_local_cell_id(int globalId)   {return globalId - this->offset;};

        Element* get_cell(int id) {return &this->elements[id];};
        pair<int, int>& get_link_face_to_cell(int id) {return this->link_face_to_cell[id];};
        pair<int, int>& get_link_face_to_vertex(int id) {return this->faces[id];};
        Node* get_centroid(int id) {return &this->centroids[id];};
        Node* get_node(int id) {return &this->nodes[id];};
        Node* get_middle_face(int id) {return &this->faceMiddlePoints[id];};
        tuple<double, double>& get_normal(int id) {return this->normals[id];};
        vector<Node>& get_nodes() {return this->nodes;};
        vector<Node>& get_centroids()  {return this->centroids;};
};      

#endif