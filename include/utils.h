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

inline bool step_1_check(double x, double y){
    if(y == 1.0)
        return true;
    return false;
}

inline bool step_2_check(double x, double y){
    if(x == 2.0)
        return true;
    return false;
}

inline bool top_check(double x, double y){
    if(y == 2) // config: backward
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
    if(x == 10) // config: backward
        return true;
    return false;
}

inline double fzero(double x, double y){
    return 0.0;
}

inline double fone(double x, double y){
    return 1.0;
}