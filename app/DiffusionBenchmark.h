double down(double x, double y) {
    return 0.0; // DIRICHLET
}

double right(double x, double y) {
    return 0.0; // DIRICHLET
}

double top(double x, double y) {
    return  1.0; // DIRICHLET
}

double left(double x, double y) {
    return 0.0; // DIRICHLET
}

double gamma(double x, double y){
    return 1.0;
}

// + Termo fonte
double source(double x, double y){
    return 0;
}
