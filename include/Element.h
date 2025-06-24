#ifndef ELEMENT_H
#define ELEMENT_H

#include "pch.h"
#include "Node.h"
#include "Edge.h"

/*
    - Tipo do elemento:
        1)  Linha
        2)  Triângulo
        3)  Quadrilatero
        15) Ponto
*/ 

class Element{
    protected:
        vector<Node*> nodes;
        vector<Edge*> edges;
    public:
        const int id;
        const int elementType; 
        
        Element(int id, int elementType, vector<Node*> nodes): id(id), elementType(elementType){
            this->nodes = nodes;
        };
        virtual ~Element() = default;

        void insert_edge(Edge* edge)            {this->edges.push_back(edge);};
        vector<Node*>& get_nodes()              {return this->nodes;};
        vector<Edge*>& get_edges()              {return this->edges;};
};

#endif