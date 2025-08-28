#include "pch.h"

double down(double x, double y) {
    return 0.0;
}

double right(double x, double y) {
    return 1.0;
}

double top(double x, double y) {
    return  0.0;
}

double left(double x, double y) {
    return 0.0;
}

double gamma(double x, double y){
    return 0.01;
}

// + Termo fonte
double source(double x, double y){
    return 0;
}

double rho(double x, double y){
    return 1.0;
}

pair<double,double> velocity(double x, double y){
    return make_pair(
       - sin(M_PI * x) * cos(M_PI * y),
       cos(M_PI * x) * sin(M_PI * y)
    );
}

// + Solução exata
double exact(double x, double y){
    return sin(M_PI * x) * sin(M_PI * y);
}