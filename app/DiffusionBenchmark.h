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
    return 1.0;
}

// + Termo fonte
double source(double x, double y){
    return -200 * x * (1-x) - 200 * y * (1-y);
}

// + Solução exata
double exact(double x, double y){
    return 100 * x * (1-x) * y * (1 - y);
}