#include "pch.h"

double down(double x, double y) {
    return 0.0; // DIRICHLET
}

double right(double x, double y) {
    return 0.0; // DIRICHLET
}

double top(double x, double y) {
    return  0.0; // DIRICHLET 
}

double left(double x, double y) {
    return 0.0; // DIRICHLET
}

double gamma(double x, double y){
    return 0.1; // difusão
}

double rho(double x, double y){
    return 1.0; // convecção
}

pair<double,double> U(double x, double y){
    return make_pair(1.0, 1.0); // velocidade
}

// + Termo fonte
double source(double x, double y){
    double pi = M_PI; 
    double term1 = pi * (std::cos(pi * x) * std::sin(pi * y) + std::sin(pi * x) * std::cos(pi * y));
    double term2 = 0.2 * pi * pi * std::sin(pi * x) * std::sin(pi * y);
    return -(term1 + term2);
}

// + Solução exata
double exact(double x, double y){
    return sin(M_PI * x) * sin(M_PI * y);
}