#include "../include/FVMSolver.h"

FVMSolver::FVMSolver(Mesh* mesh, BoundaryCondition *down, BoundaryCondition* right, BoundaryCondition* top, BoundaryCondition *left, double gamma){
    this->mesh = mesh;
    this->boundaries.push_back(down);
    this->boundaries.push_back(right);
    this->boundaries.push_back(top);
    this->boundaries.push_back(left);

    this->gamma = gamma;

    int N = this->mesh->getNumCells();
    A = new double*[N];
    _A = new double[N*N]; //declara vetor
    b = new double[N];

    /*inicialização com zeros de _A e b*/
    for(int i = 0; i < N; i++){
        b[i] = 0;
        for(int j = 0; j < N; j++){
            _A[i*N + j] = 0;
        }
    }
    for(int i = 0; i < N; i++){
        A[i] = &_A[i*N];
        b[i] = 0;
    }

    // this->printA();
}

FVMSolver::~FVMSolver(){

}

void FVMSolver::printA(){
    int N = this->mesh->getNumCells();
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            cout << A[i][j] << "\t"; 
        }
        cout << endl;
    }
}

void FVMSolver::printB(){
    
}

void FVMSolver::applyBoundariesConditions(){
}

void FVMSolver::computeA(){
    int offset = this->mesh->getTotalElements() - this->mesh->getNumCells();
    for(int i = 0; i < this->mesh->getNumCells(); i++){
        int globalCellID = this->mesh->getGlobalCellId(i);
        // Ao(icell) = 0 : já é por padrão inicializado com zero.
        // Anb(icell, :) = 0 : 0 : já são por padrão inicializados com zero.
        vector<int>* faceIdsFromCell = this->mesh->getCell(globalCellID)->getFaceIds();
        for(int j = 0; j < (*faceIdsFromCell).size(); j++){
            int gface = (*faceIdsFromCell)[j]; // global index da face
            if(this->mesh->get_link_face_to_bface(gface) != -1) // face de contorno
                continue;
            pair<int, int>* idCellsThatShareThatFace = this->mesh->get_link_face_to_cell(gface);
            int ic1 = idCellsThatShareThatFace->first; 
            int ic2 = idCellsThatShareThatFace->second;
            
            int nb = ic1 == globalCellID ? ic2 : ic1;

            A[i][globalCellID - offset] += this->gamma * this->mesh->get_face_area(gface) / this->mesh->get_deltaf(gface); // diagonal
            A[i][nb - offset] = -this->gamma *  this->mesh->get_face_area(gface) / this->mesh->get_deltaf(gface); // off-diagonal
        }
    }
    this->printA();
}

void FVMSolver::computeb(){

}