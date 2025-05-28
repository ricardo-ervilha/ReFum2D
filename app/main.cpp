#include "../include/Mesh.h"

int main(void){

    Mesh* mesh = new Mesh();
    mesh->readMesh("../inputs/malhaComplexa.msh");
    mesh->meshSummary();
    // map<int, Element>* elements = mesh->getElements();
    // cout << "#####################################################################################" << endl;
    // for(auto it = (*elements).begin(); it != (*elements).end(); ++it){
    //     if(it->second.getElementType() == 2)
    //     {
    //         cout << "Faces do elemento " << it->first << ":  ";
    //         it->second.printFaces();
    //     }    
    // }
    return 0;
}