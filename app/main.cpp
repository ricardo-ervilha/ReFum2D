#include "Mesh.h"

int main(void){
    Mesh* m = new Mesh();
    m->read_mesh("../inputs/q3x3.msh");
    return 0;
}