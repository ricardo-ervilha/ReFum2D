#include "../include/Mesh.h"

int main(void){

    Mesh* mesh = new Mesh();
    mesh->readMesh("../inputs/simpleMesh.msh");
    mesh->meshSummary();
    return 0;
}