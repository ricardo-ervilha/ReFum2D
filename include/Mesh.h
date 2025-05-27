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

using namespace std;

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

        /*Faces*/
        int nfaces; //Qtd total de faces
        vector<pair<int, int>> faces; //Vetor de faces. Ordem das faces já vem da indexação do vetor.
        // padrão (min, max), isto é: face 1,2 = face 2,1, mas no vetor só tem 1,2.

        /*Geral*/
        int geom_type; //Tipo da geometria (2d ou 3d)
        int nbfaces; //Qtd total de faces de contorno


        map<int, Node> centroids; //int nodeID => centroide da célula
        vector<double> faceAreas; // como está tudo 2D por enquanto, aqui terá o comprimento somente.
        vector<Node> faceMiddlePoints; //armazena o ponto central das faces. 

        map<int,int> elementTypeToNumNodes; // int id do elemento => int número de nós que o compõem

        /*Funções privadas*/
        void preProcessing();
    public:
        Mesh();
        ~Mesh();
        void readMesh(string filepath); // Função para ler os dados da malha .msh
        void meshSummary();
        int getGeomType() {return this->geom_type;};
        int getNFaces()  {return this->nfaces;};
        int getNbFaces() const { return nbfaces; };
        int getNNodes() const { return nnodes; };
        int getTotalElements() const { return totalElements; };
        int getTotalPhysicalEntities() const { return totalPhysicalGroups; };
        vector<Node>* getNodes() { return &nodes; };
        map<int, PhysicalGroup>* getPhysicalGroups() { return &physicalGroups; };
        vector<Element>* getElements() { return &elements; };
        map<int, Node>* getCentroids() { return &centroids; };
        vector<pair<int, int>>* getFaces() { return &faces; };
        vector<double>* getFaceAreas() { return &faceAreas; };
        vector<Node>* getFaceMiddlePoints() { return &faceMiddlePoints; };
        map<int, int>* getElementTypeToNumNodes() { return &elementTypeToNumNodes; };
};

#endif