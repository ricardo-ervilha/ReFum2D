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

/*Função hash para o unordered set*/
struct CantorPairingFunction {
    // injetora para inteiros >= 0.
    std::size_t operator()(const std::pair<int, int>& p) const noexcept {
        int a = p.first;
        int b = p.second;
        return (a + b) * (a + b + 1) / 2 + b;
    }
};

struct Edge {
    int from, to; // mantém a ordem original

    Edge(int a, int b) : from(a), to(b) {}

    pair<int, int> asKey() const {
        return (from < to) ? make_pair(from, to) : make_pair(to, from);
    }
};

class Mesh{
    private:
        /*Nodes*/
        vector<Node> nodes; //guarda os nós 
        int nnodes; //Qtd total de nós/vértices

        /*Elements*/
        vector<Element> elements; //guarda os elementos (células e outros que vem junto)
        int totalElements; //Qtd total de elementos 

        /*Grupos físicos*/
        map<int, PhysicalGroup> physicalGroups; //mapeamento: int PhysicalGroupID => Grupo Físico
        int totalPhysicalGroups; //Qtd total de grupos físicos

        /*Células*/
        int ncells; //Qtd total de células
        vector<int> cells; //guarda o id das células.
        vector<double> volumes; //volumes de cada célula

        /*Faces*/
        int nfaces; //Qtd total de faces
        vector<pair<int, int>> faces; //Vetor de faces. Ordem das faces já vem da indexação do vetor. | padrão (min, max), isto é: face 1,2 = face 2,1, mas no vetor só tem 1,2.
        vector<pair<int,int>> link_face_to_cell; // para a i-ésima face, o pair contém as células que compartilham a face.

        /*Geral*/
        int geom_type; //Tipo da geometria (2d ou 3d)
        map<int,int> elementTypeToNumNodes; // int id do elemento => int número de nós que o compõem
        
        /*Faces de contorno*/
        int nbfaces; //Qtd total de faces de contorno
        vector<int> link_face_to_bface; //vetor que acompanha a indexação das faces, e em cada entrada correspondendo a id face global, diz qual id da bface, ou retorna -1 caso não seja boundary face.
        vector<int> link_bface_to_face; // vetor do tamanho de faces de contorno. Para cada indice do vetor, diz qual o indice global da boundary face.
        unordered_map<pair<int,int>, int, CantorPairingFunction> facesPairToKey; // dado o par de chaves retorna o índice => vai ser útil na hora das faces de contorno

        /*Centroides*/
        vector<Node> centroids; // acompanha a indexação das células

        /*Normais*/
        vector<tuple<double, double>> normals; // normais de cada face. Acompanha a indexação de faces
        
        /*Área da Face (3d) ou Comprimento da face (1D)*/
        vector<double> faceAreas; // acompanha a indexação do faces.
        
        /*Ponto médio das faces*/
        vector<Node> faceMiddlePoints; //armazena o ponto central das faces. 

        /*Funções privadas*/
        void preProcessing(); //Calcula valores de diversas coisas que não foram feitos no readMesh

        /*Deltaf*/
        vector<double> deltafs; // acompanha a indexação da face (distancia entre os centroids das celulas que compartilham a face)
    public:
        Mesh();
        ~Mesh();
        void readMesh(string filepath); // Função para ler os dados da malha .msh
        void meshSummary();
        void printInfoGeral();
        void printInfoNodes();
        void printInfoElements();
        void printInfoPhyisicalGroups();
        void printInfoCells();
        void printInfoFaces();
        void printInfoCentroids();
        void printInfoLinkFaceToCell();
        void printinfoNormalSigns();
        void printInfoNormalValues();
        void printInfoVolumes();
        void printInfoLinkFaceToBface();
        void printInfoLinkBfaceToFace();
        void printDeltaFs();

        int getNumCells() {return this->ncells;};
        int getTotalElements()   {return this->totalElements;};
        int getGlobalCellId(int i) {return this->cells[i];};
        Element* getCell(int id) {return &this->elements[id];};
        int get_link_face_to_bface(int id) {return this->link_face_to_bface[id];};
        pair<int, int>* get_link_face_to_cell(int id) {return &this->link_face_to_cell[id];};
        double get_face_area(int id) {return this->faceAreas[id];};
        double get_deltaf(int id) {return this->deltafs[id];};
        pair<int, int>* get_id_of_nodes_that_make_face(int id) {return &this->faces[id];};
        Node* get_centroid(int id) {return &this->centroids[id];};
        Node* get_node(int id) {return &this->nodes[id];};
};

#endif