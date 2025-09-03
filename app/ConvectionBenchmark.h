#include "pch.h"
double down(double x, double y) {
    return 0.0;
}

double right(double x, double y) {
    return 0.0;
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
    return - (0.02 * M_PI * M_PI * sin(M_PI * x) * sin(M_PI * y) + \
            M_PI * (cos(M_PI * x) * sin(M_PI * y) + 0.5 * sin(M_PI * x) * cos(M_PI * y)));
}

double rho(double x, double y){
    return 1.0;
}

pair<double,double> velocity(double x, double y){
    return make_pair(
        1,
        0.5
    );
}

// + Solução exata
double exact(double x, double y){
    return sin(M_PI * x) * sin(M_PI * y);
}