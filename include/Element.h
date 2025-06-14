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
        vector<int> nodes; //id dos nós que compõem o Elemento.
        vector<int> facesIds; //id das faces que compõem o Elemento.
        vector<int> normalSigns; 
        /*
            Sobre o normalSigns:
                1: já está apontando pra fora do elemento ou -1 para "consertar" e apontar para fora do elemento.
        */

    public:
        Element(int elementType, vector<int> nodes) {
            this->elementType = elementType;
            this->nodes = nodes;
        };
        ~Element() {};
        
        int get_element_type()    {return this->elementType;};
        vector<int>& get_normal_sign() {return this->normalSigns;};
        vector<int>& get_nodes() {return this->nodes;};
        vector<int>& get_face_ids() {return this->facesIds;};
        
        void insert_face(int idFace) {this->facesIds.push_back(idFace);};
        void insert_normal_sign(int normalSign) {this->normalSigns.push_back(normalSign);};
};

#endif