#pragma once

#include "pch.h"

/*Função para computar a distância entre dois pontos.*/
inline double distance(double x1, double y1, double x2, double y2) {
    double dx = x1 - x2;
    double dy = y1 - y2;
    return sqrt(dx * dx + dy * dy);
}