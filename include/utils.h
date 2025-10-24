#pragma once

#include "pch.h"

/*Função para computar a distância entre dois pontos.*/
inline double distance(double x1, double y1, double x2, double y2) {
    double dx = x1 - x2;
    double dy = y1 - y2;
    return sqrt(dx * dx + dy * dy);
}

inline int get_neighbor(pair<int,int>& id_cells_share_face, int P){
    int ic1 = id_cells_share_face.first; 
    int ic2 = id_cells_share_face.second;

    return ic1 == P ? ic2 : ic1;
}

inline bool circle_check(double x, double y){
    // Se o ponto (x,y) pertence ao círculo então: sqrt((x-xc)**2 + (y-yc)**2) = r ou sqrt((x-xc)**2 + (y-yc)**2) - r = 0
    // como o gmsh gera refinando o polígono, nao sera exatamente zero, ai estou pedindo pra ser menor q uma tolerância.
    // círculo está centrado em 0.2 e 0.2
    // raio dele é 0.05
    if(fabs(sqrt(pow(x-0.2,2) + pow(y-0.2,2))- 0.05) <= 1e-3)
        return true;
    else
        return false;
}

inline bool step_1_check(double x, double y){
    if(y == 0.5)
        return true;
    return false;
}

inline bool step_2_check(double x, double y){
    if(x == 0.5)
        return true;
    return false;
}

inline bool top_check(double x, double y){
    // if(y == 1) // config: backward
    // if(y == 0.41) // config: flow_over_cylinder
    if(y == 1) // config: lid
        return true;
    return false;
}

inline bool bottom_check(double x, double y){
    if(y == 0.0)
        return true;
    return false;
}

inline bool left_check(double x, double y){
    if(x == 0.0)
        return true;
    return false;
}

inline bool right_check(double x, double y){
    // if(x == 10) // config: backward
    // if(x == 2.2) // config: flow over cylinder
    if(x == 1) // config: lid
        return true;
    return false;
}

inline double fzero(double x, double y){
    return 0.0;
}

inline double fone(double x, double y){
    return 1.0;
}