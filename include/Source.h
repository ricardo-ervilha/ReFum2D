#ifndef SOURCE_H
#define SOURCE_H
#include "FVMSolver.h"

class Source{
    private:
        // & Função que recebe x, y para pegar o valor do termo fonte naquele vértice da malha.
        double (*sourcefunc)(double, double); 

        // # guarda o solver para recuperar malha e ed's
        FVMSolver* solver;

        // & Interpola os termos fontes em sourcecc.
        void interpolate_sources();
    public:
         // # Irá armazenar o valor de do termo fonte em cada centro da celula 
        arma::vec sourcecc;

        Source(FVMSolver* s, double (*source)(double, double));
        ~Source();

        void assemblyCoefficients();
};

#endif