#ifndef NODE_H
#define NODE_H

#include <vector>

using namespace std;
class Node{
    private:
        /*coordenadas*/
        double x; 
        double y;
        double z;
        vector<double> distanceCentroids; //guarda a distância para o centroide
        vector<int> idCellRelativeToCentroid; //guarda para qual célula correspondente é aquela distância.
        /*
            exemplo de possíveis vetores:
            - 0.25, 0.3, 0.4
            - 15, 19, 24
            Assim é possível ir iterando para preencher a matriz futuramente.
        */
    public:
        Node(double x, double y, double z)  
        {
            this->x = x; 
            this-> y = y; 
            this->z = z;
        };
        ~Node() {};
        double getX() {return this->x;};
        double getY() {return this->y;};
        double getZ() {return this->z;};
        void insertDistanceCentroids(double distance) {this->distanceCentroids.push_back(distance);};
        void insertIdCellRelativeToCentroid(int id) {this->idCellRelativeToCentroid.push_back(id);};
        vector<double>* getDistanceCentroids() {return &this->distanceCentroids;};
        vector<int>* getIdCellRelativeToCentroid() {return &this->idCellRelativeToCentroid;};
};

#endif

