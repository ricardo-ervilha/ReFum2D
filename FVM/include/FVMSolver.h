#ifndef FVMSOLVER_H
#define FVMSOLVER_H

#include <iostream>
#include <vector>
#include "Node.h"

using namespace std;

class FVMSolver{
    private:
        int geom_type; //Tipo da geometria (2d ou 3d)
        int ncells; //Qtd total de células
        int nfaces; //Qtd total de faces
        int nbfaces; //Qtd total de faces de contorno
        int nnodes; //Qtd total de nós/vértices

        vector<Node> nodes; //Vetor de nodes, contendo suas posições em x, y e z.
    public:
        FVMSolver();
        ~FVMSolver();
        void readMesh(string filepath); // Função para ler os dados da malha .msh
};

#endif