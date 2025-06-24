#ifndef NODE_H
#define NODE_H

#include "pch.h"

using namespace std;
class Node{
    public:
        const int id;
        const double x; 
        const double y;
        Node(int id, double x, double y): id(id), x(x), y(y) {};
        ~Node() {};
};

#endif

