#ifndef NODE_H
#define NODE_H

class Node{
    private:
        double x;
        double y;
        double z;
    public:
        Node(double x, double y, double z)  {this->x = x; this-> y = y; this->z = z;};
        ~Node() {};
        double getX() {return this->x;};
        double getY() {return this->y;};
        double getZ() {return this->z;};
};

#endif

