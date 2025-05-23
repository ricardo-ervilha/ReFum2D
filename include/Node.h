#ifndef NODE_H
#define NODE_H

class Node{
    private:
        int id;
        double x;
        double y;
        double z;
    public:
        Node(int id, double x, double y, double z)  
        {
            this->id = id;
            this->x = x; 
            this-> y = y; 
            this->z = z;
        };
        Node(){
            this->id = 0;
            this->x = 0;
            this->y = 0;
            this->z = 0;
        }
        ~Node() {};
        int getId()   {return this->id;};
        double getX() {return this->x;};
        double getY() {return this->y;};
        double getZ() {return this->z;};
};

#endif

