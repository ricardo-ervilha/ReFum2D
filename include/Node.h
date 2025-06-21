#ifndef NODE_H
#define NODE_H

#include "pch.h"

using namespace std;
class Node{
    private:
        double x; 
        double y;
        double z;

        vector<double> distance_centroids; //guarda a distância para o centroide
        vector<int> id_cell_relative_to_centroids; //guarda para qual célula correspondente é aquela distância.
        /*
            #Exemplo de possíveis vetores:
            -> 0.25, 0.3, 0.4
            -> 15, 19, 24
            #Significam que esse nó faz parte das células 15, 15 e 24; e as distâncias desse nó aos centróides dessas células respectivamente são 0.25, 0.3 e 0.4.
        */
    public:
        Node(double x, double y, double z)  
        {
            this->x = x; 
            this-> y = y; 
            this->z = z;
        };
        ~Node() {};
        
        double get_x() {return this->x;};
        double get_y() {return this->y;};
        double get_z() {return this->z;};

        vector<double>& get_distance_centroids() {return this->distance_centroids;};
        vector<int>& get_id_relative_to_centroids() {return this->id_cell_relative_to_centroids;};
        
        void insert_distance_centroids(double distance) {this->distance_centroids.push_back(distance);};
        void insert_id_relative_to_centroid(int id) {this->id_cell_relative_to_centroids.push_back(id);};
};

#endif

