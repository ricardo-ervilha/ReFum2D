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
    skew = new double[N];
    u = new double[N];

    /*inicialização com zeros de _A e b*/
    for(int i = 0; i < N; i++){
        b[i] = 0;
        skew[i] = 0;
        u[i] = 0; //inicializa a solução com 0
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
    cout << "[" << endl;
    for (int i = 0; i < N; i++)
    {
        cout << " [";
        for (int j = 0; j < N; j++)
        {
            cout << setw(8) << fixed << setprecision(2) << A[i][j];
            if (j < N - 1)
                cout << ", ";
        }
        cout << "]" << (i < N - 1 ? "," : "") << endl;
    }
    cout << "]" << endl;
}

void FVMSolver::printB(){
    int N = this->mesh->getNumCells();
    cout << " [";
    for(int i = 0; i < N; i++){
        cout << setw(8) << fixed << setprecision(2) << b[i];
    }
    cout << "]" << endl;
    cout << endl;
    cout << endl;
}

void FVMSolver::applyBoundariesConditions(){
    int nbfaces = this->mesh->get_num_boundary_faces();
    int offset = this->mesh->getTotalElements() - this->mesh->getNumCells();
    for(int i = 0; i < nbfaces; i++){
        int gface = this->mesh->get_link_bface_to_face(i);
        int ic = this->mesh->get_link_face_to_cell(gface)->first;
        b[ic - offset] += this->gamma * this->mesh->get_face_area(gface) * 10 / this->mesh->get_deltaf(gface);
    }
    this->printB();
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

double FVMSolver::phiv(Node* n){
    int offset = this->mesh->getTotalElements() - this->mesh->getNumCells(); //offset para dar id local
    vector<double>* distanceCentroids = n->getDistanceCentroids();
    vector<int>* idCellRelativeToCentroid =  n->getIdCellRelativeToCentroid(); //retorna id global
    double phiv = 0;
    for(int i = 0; i < distanceCentroids->size(); i++){
        phiv += (*distanceCentroids)[i] * this->u[(*idCellRelativeToCentroid)[i] - offset];        
    }    
    return phiv;
}

void FVMSolver::computeb(){
    int offset = this->mesh->getTotalElements() - this->mesh->getNumCells();
    for(int i = 0; i < this->mesh->getNumCells(); i++){
        int globalCellID = this->mesh->getGlobalCellId(i);
        vector<int>* faceIdsFromCell = this->mesh->getCell(globalCellID)->getFaceIds(); // id das faces daquela célula
        for(int j = 0; j < (*faceIdsFromCell).size(); j++){
            int gface = (*faceIdsFromCell)[j]; // global index da face
            
            if(this->mesh->get_link_face_to_bface(gface) != -1) // face de contorno
                continue;
            
            //encontrando os indices dos vertices que compõem essa face
            pair<int, int>* idNodes = this->mesh->get_id_of_nodes_that_make_face(gface);
            int va = idNodes->first;
            int vb = idNodes->second;

            pair<int, int>* idCellsThatShareThatFace = this->mesh->get_link_face_to_cell(gface);
            int ic1 = idCellsThatShareThatFace->first; 
            int ic2 = idCellsThatShareThatFace->second;

            Node* ic1El = this->mesh->get_centroid(ic1 - offset);
            Node* ic2El = this->mesh->get_centroid(ic2 - offset);

            double Lf1 = ic2El->getX() - ic1El->getX();
            double Lf2 = ic2El->getY() - ic2El->getY();

            Node* vaN = this->mesh->get_node(va);
            Node* vbN = this->mesh->get_node(vb);

            double tf1 = vaN->getX() - vbN->getX();
            double tf2 = vaN->getY() - vbN->getY();
            
            // tangente unitária
            tf1 = tf1 / this->mesh->get_face_area(gface);
            tf2 = tf2 / this->mesh->get_face_area(gface);
            
            double utf_dot_Lf = tf1 * Lf1 + tf2 * Lf2;

            skew[globalCellID - offset] += -this->gamma * utf_dot_Lf * (phiv(vaN) - phiv(vbN)) / this->mesh->get_deltaf(gface);
        }
    }
    this->printB();
}