#ifndef ELEMENT_H
#define ELEMENT_H

#include "Node.h"
#include <vector>

using namespace std;

class Element{
    private:
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
        vector<int> facesIds; //id das faces que compõem o element.
        vector<int> normalSigns; // 1: já está apontando pra fora | -1 conserta para apontar pra fora
    public:
        Element(int elementType, vector<int> nodes) {
            this->elementType = elementType;
            this->nodes = nodes;
        };
        ~Element() {};
        int getElementType()    {return this->elementType;};
        vector<int>* getNodes() {return &this->nodes;};
        void insertFace(int idFace) {this->facesIds.push_back(idFace);};
        vector<int>* getFaceIds() {return &this->facesIds;};
        void insertNormalSign(int normalSign) {this->normalSigns.push_back(normalSign);};
        vector<int>* getNormalSigns() {return &this->normalSigns;};
};

#endif