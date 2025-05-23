#ifndef MESH_H
#define MESH_H

#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <tuple>
#include <math.h>
#include "Node.h"
#include "PhysicalGroup.h"
#include "Element.h"

using namespace std;

class Mesh{
    private:
        /*Variáveis e atributos*/
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

        map<int, Node> centroids; //int nodeID => centroide da célula
        vector<pair<int, int>> faces; // padrão (min, max), isto é: face 1,2 = face 2,1, mas no vetor só tem 1,2.
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
        int getTotalPhysicalEntities() const { return totalPhysicalEntities; };
        map<int, Node>* getNodes() { return &nodes; };
        map<int, PhysicalGroup>* getPhysicalGroups() { return &physicalGroups; };
        map<int, Element>* getElements() { return &elements; };
        map<int, Node>* getCentroids() { return &centroids; };
        vector<pair<int, int>>* getFaces() { return &faces; };
        vector<double>* getFaceAreas() { return &faceAreas; };
        vector<Node>* getFaceMiddlePoints() { return &faceMiddlePoints; };
        map<int, int>* getElementTypeToNumNodes() { return &elementTypeToNumNodes; };
};

#endif