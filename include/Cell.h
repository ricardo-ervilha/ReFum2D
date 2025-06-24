#ifndef CELL_H
#define CELL_H

#include "pch.h"
#include "Element.h"
#include "Edge.h"

class Cell : public Element{
    private:
        vector<int> nsigns;
        double area;
        pair<double, double> centroid;
        void compute_properties(){
            // centroide
            double xc = 0, yc = 0;
            for(int i = 0; i < this->nodes.size(); i++){
                xc += nodes[i]->x;
                yc += nodes[i]->y;
            }
            xc = xc/this->nodes.size(); yc = yc/this->nodes.size(); 
            centroid = make_pair(xc,yc);
        };
    public: 
        const int idCell;
        Cell(int id, int elementType, vector<Node*> nodes, int idCell) : Element(id, elementType, nodes), idCell(idCell) {
            this->compute_properties();
        };
        ~Cell()     {};

        void insert_nsign(int nsign) {this->nsigns.push_back(nsign);};
        void set_area(double area)  {this->area = area;};
        
        pair<double,double>& get_centroid() {return this->centroid;};

        vector<int>& get_nsigns()   {return this->nsigns;};
        double get_area()   {return this->area;};
};

#endif