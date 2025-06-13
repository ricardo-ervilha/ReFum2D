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
    Qs = new double[N];

    /*inicialização com zeros de _A e b*/
    for(int i = 0; i < N; i++){
        b[i] = 0;
        skew[i] = 0;
        u[i] = 0; //inicializa a solução com 0
        Qs[i] = 0;
        for(int j = 0; j < N; j++){
            _A[i*N + j] = 0;
        }
    }
    for(int i = 0; i < N; i++){
        A[i] = &_A[i*N];
        b[i] = 0;
    }

    this->computeQs();
}

FVMSolver::~FVMSolver(){

}

double FVMSolver::calculateQ(double x, double y){
    return -200*x*(1-x) - 200*y*(1-y);
}

void FVMSolver::computeQs(){
    for(int i = 0; i < this->mesh->getNumCells(); i++){
        int globalCellId = this->mesh->getGlobalCellId(i);
        Element *cell = this->mesh->getCell(globalCellId);
        vector<int>* nodeIds = cell->getNodes();

        Node* n1 = this->mesh->get_node((*nodeIds)[0]);
        Node* n2 = this->mesh->get_node((*nodeIds)[1]);
        Node* n3 = this->mesh->get_node((*nodeIds)[2]);

        double Q1 = this->calculateQ(n1->getX(), n1->getY());
        double Q2 = this->calculateQ(n2->getX(), n2->getY());
        double Q3 = this->calculateQ(n3->getX(), n3->getY());

        this->Qs[i] = (Q1 + Q2 + Q3)/3.0;
    }
}

void FVMSolver::printA(){
    int N = this->mesh->getNumCells();
    cout << "[" << endl;
    for (int i = 0; i < N; i++)
    {
        cout << " [";
        for (int j = 0; j < N; j++)
        {
            cout << setw(8) << fixed << setprecision(3) << A[i][j];
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
        cout << setw(8) << fixed << setprecision(3) << b[i];
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
            {   
                Node* centroid = this->mesh->get_centroid(globalCellID - offset);

                Node* middleOfTheFace = this->mesh->get_middle_points(gface);

                double lb1 = middleOfTheFace->getX() - centroid->getX();
                double lb2 = middleOfTheFace->getY() - centroid->getY();

                tuple<double, double>* normals = this->mesh->get_normal(gface);

                double deltaB = lb1 * get<0>(*normals) + lb2 * get<1>(*normals);

                A[i][globalCellID - offset] += this->gamma * this->mesh->get_face_area(gface) / deltaB;

            }else{
                pair<int, int>* idCellsThatShareThatFace = this->mesh->get_link_face_to_cell(gface);
                int ic1 = idCellsThatShareThatFace->first; 
                int ic2 = idCellsThatShareThatFace->second;
                
                int nb = ic1 == globalCellID ? ic2 : ic1;

                A[i][globalCellID - offset] += this->gamma * this->mesh->get_face_area(gface) / this->mesh->get_deltaf(gface); // diagonal
                A[i][nb - offset] = -this->gamma *  this->mesh->get_face_area(gface) / this->mesh->get_deltaf(gface); // off-diagonal
            }
        }
    }
    // this->printA();
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
        b[globalCellID - offset] = skew[globalCellID - offset] - Qs[globalCellID - offset]*this->mesh->get_volume(globalCellID - offset);
    }
    this->printB();
}

double* FVMSolver::copyVector(double* M, int N){
    double* X = new double[N];

    for(int i = 0; i < N; i++)  
        X[i] = M[i];

    return X;
}

double FVMSolver::max_norm_difference(double* u, double* u_old, int N){
    /*Computa a norma do máximo da diferença entre dois vetores*/
    double max_diff = 0.0;
    for (int i = 0; i < N; ++i) {
        double diff = std::abs(u[i] - u_old[i]);
        if (diff > max_diff) {
            max_diff = diff;
        }
    }
    return max_diff;
}

void FVMSolver::GaussSeidel(double tol) {
    double error = 1.0;
    int N = this->mesh->getNumCells();
    double sum;
    int iter = 0; 
    // Cria uma cópia da solução inicial para cálculo do erro
    double* tmp = this->copyVector(this->u, N);

    while (error > tol) {
        // Salva a solução anterior em tmp
        for (int i = 0; i < N; i++) {
            tmp[i] = u[i];
        }

        for (int i = 0; i < N; i++) {
            sum = 0.0;
            for (int j = 0; j < N; j++) {
                if (j != i)
                    sum += A[i][j] * u[j];
            }
            u[i] = (b[i] - sum) / A[i][i];
        }
        
        iter += 1;
        error = max_norm_difference(u, tmp, N);
        std::cout << "Valor do erro: " << error << std::endl;
    }

    std::cout << "Solução final: " << std::endl;
    std::cout << " [";
    for (int i = 0; i < N; i++) {
        std::cout << std::setw(8) << std::fixed << std::setprecision(2) << u[i];
    }
    std::cout << "]" << std::endl << std::endl;

    delete[] tmp;  // Libera a memória alocada para tmp
}

void FVMSolver::saveSolution(){
    string filename = "../inputs/solution.vtk";
    ofstream vtk_file(filename);

    if (!vtk_file.is_open()) {
        std::cerr << "Erro ao criar o arquivo '" << filename << "'.\n";
        exit(1);
    }

    vtk_file << "# vtk DataFile Version 3.0\n";
    vtk_file << "Malha 2D com triângulos e temperatura\n";
    vtk_file << "ASCII\n";
    vtk_file << "DATASET UNSTRUCTURED_GRID\n\n";

    vtk_file << "POINTS " << this->mesh->getNumCells()  << " double" << endl;
    vector<Node>* nodes = this->mesh->getNodes();
    for(int i = 0; i < this->mesh->getNumCells(); i++){
        vtk_file << (*nodes)[i].getX() << " " << (*nodes)[i].getY() << " " << (*nodes)[i].getZ() << endl;
    }
    vtk_file << endl;

    vtk_file << "CELLS " << this->mesh->getNumCells() << " " <<  4 * this->mesh->getNumCells() << endl;
    for(int i = 0; i < this->mesh->getNumCells(); i++){
        int gcell = this->mesh->getGlobalCellId(i);
        Element* e = this->mesh->getCell(gcell);
        vector<int>* nodeIds = e->getNodes();
        vtk_file << "3" << " " << (*nodeIds)[0] << " " << (*nodeIds)[1] << " " << (*nodeIds)[2] << endl;
    }

    vtk_file << endl;

    vtk_file << "CELL_TYPES " << this->mesh->getNumCells() << endl;
    for(int i = 0; i < this->mesh->getNumCells(); i++){
        vtk_file << "5" << endl;
    }

    vtk_file << endl;
    vtk_file << "CELL_DATA " << this->mesh->getNumCells() << endl;
    vtk_file << "SCALARS Temperatura double 1\n";
    vtk_file << "LOOKUP_TABLE default\n";
    for(int i = 0; i < this->mesh->getNumCells(); i++){
        vtk_file << this->u[i] << endl;
    }

    vtk_file.close();
}

void FVMSolver::computeErrorExact(){
    vector<Node>* centroids = this->mesh->getCentroids();
    int n = this->mesh->getNumCells();
    double* exact = new double[n]; 
    for(int i = 0; i < (*centroids).size(); i++){
        double x = (*centroids)[i].getX();
        double y = (*centroids)[i].getY();
        exact[i] = 100*x*(1-x)*y*(1-y);
    }
    double error = this->max_norm_difference(exact, this->u, n);
    cout << "Erro máximo: " << error << endl;
}