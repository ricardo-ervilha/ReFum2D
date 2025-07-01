#ifndef EDGE_H
#define EDGE_H

#include "pch.h"
#include "utils.h"
#include "Node.h"
#include <limits>

class Edge{
    private:
        double length;
        double df; /*distância entre os centros das celulas que dividem a aresta.*/
        pair<double,double> middle;
        pair<double,double> normal;
        pair<int,int> link_face_to_cell;

        void compute_properties(){
            length = distance(from->x, from->y, to->x, to->y);
            middle = make_pair((from->x + to->x)*0.5, (from->y + to->y)*0.5);
            normal = make_pair((to->y - from->y)/length, (from->x - to->x)/length);
        };
    public:
        const int id;
        const Node* from;
        const Node* to;
        
        Edge(int id, Node* from, Node* to) : id(id), from(from), to(to) {
            this->compute_properties();
            this->link_face_to_cell = make_pair(-1, -1); // inicializa com flag
        };
        ~Edge() {};

        void set_link_face_to_cell(int lftc, int pos) 
        {
            if(pos == 1)
                this->link_face_to_cell.first = lftc;
            else if(pos == 2)
                this->link_face_to_cell.second = lftc;
        };
        pair<int,int>& get_link_face_to_cell()    {return link_face_to_cell;};

        double get_df()     {return this->df;};
        void set_df(double df)     {this->df = df;};

        double get_length()                   {return this->length;};
        pair<double,double>& get_middle()     {return this->middle;};
        pair<double,double>& get_normal()     {return this->normal;};

        bool is_boundary_face() {
            if(link_face_to_cell.first == -1 || link_face_to_cell.second == -1)
                return true;
            else
                return false;
        }

};

#endif