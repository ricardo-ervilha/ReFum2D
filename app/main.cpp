#include "../include/Mesh.h"

int main(void){

    Mesh* mesh = new Mesh();
    mesh->readMesh("../inputs/complexMesh.msh");
    mesh->meshSummary();
    return 0;
}