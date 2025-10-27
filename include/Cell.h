#ifndef CELL_H
#define CELL_H

#include "Edge.h"

class Cell{
    private:
        vector<int> nsigns;
        double area;
        pair<double, double> centroid;
        vector<Node*> nodes;
        vector<Edge*> edges;

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
        const int id;
        const int cellType;
        Cell(int id, int cellType, vector<Node*> nodes) : id(id), cellType(cellType) {
            this->nodes = nodes;
            this->compute_properties();
        };
        ~Cell()     {};

        void insert_nsign(int nsign) {this->nsigns.push_back(nsign);};
        void set_area(double area)  {this->area = area;};
        
        pair<double,double>& get_centroid() {return this->centroid;};

        vector<int>& get_nsigns()   {return this->nsigns;};
        double get_area()   {return this->area;};

        void insert_edge(Edge* edge)            {this->edges.push_back(edge);};
        vector<Node*>& get_nodes()              {return this->nodes;};
        vector<Edge*>& get_edges()              {return this->edges;};

};

#endif