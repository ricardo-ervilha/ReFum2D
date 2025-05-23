#ifndef ELEMENT_H
#define ELEMENT_H

#include "Node.h"
#include <vector>

using namespace std;

class Element{
    private:
        int id;
        /*
            Tipo do elemento:
            1)  Linha
            2)  Triângulo
            3)  Quadrilatero
            4)  Tetraedro
            15) Ponto
        */ 
        int elementType; 
        vector<int> nodes; //id dos nós que compõem o element.
        vector<int> faces; //id das faces que compõem o element.
    public:
        Element(int id, int elementType, vector<int> nodes) {
            this->id = id;
            this->elementType = elementType;
            this->nodes = nodes;
        };
        ~Element() {};
        int getId()             {return this->id;};
        int getElementType()    {return this->elementType;};
        vector<int>* getNodes() {return &this->nodes;};
        void insertFace(int idFace) {this->faces.push_back(idFace);};
        void printNodes(){
            for(int i = 0; i < nodes.size(); i++){
                cout << nodes[i] << " ";
            }
            cout << endl;
        };
        void printFaces(){
            for(int i = 0; i < faces.size(); i++){
                cout << faces[i] << " ";
            }
            std::cout << std::endl;
        }
};

#endif