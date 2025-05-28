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
        int nbfaces; //Qtd total de faces de contorno
        map<int,int> elementTypeToNumNodes; // int id do elemento => int número de nós que o compõem

        /*Centroides*/
        vector<Node> centroids; // acompanha a indexação das células

        /*Normais*/
        vector<tuple<double, double>> normals; // normais de cada face. Acompanha a indexação de faces
        
        /*Área da Face (3d) ou Comprimento da face (1D)*/
        vector<double> faceAreas; // acompanha a indexação do faces.
        
        /*Ponto médio das faces*/
        vector<Node> faceMiddlePoints; //armazena o ponto central das faces. 

        /*Conectividades*/
        //link_cell_to_face | Já implementado através do vetor de faces da célula.
        //link_face_to_node | Já implementado através da face guardar o pair de nós
        //link_cell_to_node | Já implementado através do vetor de nós da célula.
        //link_face_to_bface|
        //link_bface_to_face|

        /*Funções privadas*/
        void preProcessing(); //Calcula valores de diversas coisas que não foram feitos no readMesh
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
        vector<Node>* getCentroids() { return &centroids; };
        vector<pair<int, int>>* getFaces() { return &faces; };
        vector<double>* getFaceAreas() { return &faceAreas; };
        vector<Node>* getFaceMiddlePoints() { return &faceMiddlePoints; };
        map<int, int>* getElementTypeToNumNodes() { return &elementTypeToNumNodes; };
};

#endif