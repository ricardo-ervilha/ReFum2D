#include "Mesh.h"

int main(void){

    Mesh* mesh = new Mesh();
    mesh->readMesh("../inputs/malhaSimples.msh");
    mesh->meshSummary();
    return 0;
}