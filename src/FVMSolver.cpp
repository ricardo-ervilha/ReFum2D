#include "../include/FVMSolver.h"

FVMSolver::FVMSolver(Mesh* mesh, BoundaryCondition *down, BoundaryCondition* right, BoundaryCondition* top, BoundaryCondition *left, double (*g)(double, double), double (*sourceTerm)(double, double)){
    this->mesh = mesh;

    /*Salva as condições usando a mesma orientação definida.*/
    this->boundaries.push_back(down);
    this->boundaries.push_back(right);
    this->boundaries.push_back(top);
    this->boundaries.push_back(left);

    this->gammafunc = g;
    this->sourcefunc = sourceTerm;

    int ncells = this->mesh->get_num_cells();

    /*Inicialização das EDs associadas a resolução do problema...*/
    this->A = vector<vector<double>>(ncells, vector<double>(ncells, 0)); // ncells x ncells : 0's
    this->b = vector<double>(ncells, 0); // 1 x n : 0's
    this->u = vector<double>(ncells, 0); // 1 x n: 0's
    this->source = vector<double>(ncells, 0); // 1 x n: 0's
    this->gamma = vector<double>(ncells, 0); // 1 x n : 0's

    this->pre_processing();
}

FVMSolver::~FVMSolver(){
}

void FVMSolver::pre_processing(){
    int ncells = mesh->get_num_cells();
    int gcell;
    Element* cell;
    double source1, source2, source3;
    double gamma1, gamma2, gamma3;
    for(int i = 0; i < ncells; i++){ // itera pelos índices locais
        gcell = mesh->get_global_cell_id(i);
        cell = mesh->get_cell(gcell);
        vector<int>& nodeIds = cell->get_nodes();

        Node* n1 = mesh->get_node(nodeIds[0]);
        Node* n2 = mesh->get_node(nodeIds[1]);
        Node* n3 = mesh->get_node(nodeIds[2]);

        source1 = this->sourcefunc(n1->get_x(), n1->get_y());
        source2 = this->sourcefunc(n2->get_x(), n2->get_y());
        source3 = this->sourcefunc(n3->get_x(), n3->get_y());

        gamma1 = this->gammafunc(n1->get_x(), n1->get_y());
        gamma2 = this->gammafunc(n2->get_x(), n2->get_y());
        gamma3 = this->gammafunc(n3->get_x(), n3->get_y());

        /*Interpolação das fontes é só usando uma média mesmo dos valores nos nós.*/
        this->source[i] = (source1 + source2 + source3) / 3;

        /*Interpolação do gamma é só usando a média por enquanto.*/
        this->gamma[i] = (gamma1 + gamma2 + gamma3) / 3;
    }
}

double FVMSolver::phiv(Node* n){
    vector<double>& distCentroids = n->get_distance_centroids();
    vector<int>& idRelativeCentroids = n->get_id_relative_to_centroids();

    double phiv = 0;

    for(int i = 0; i < distCentroids.size(); i++){
        int lcell = mesh->get_local_cell_id(idRelativeCentroids[i]);
        phiv += distCentroids[i] * u[lcell];
    }

    return phiv;
}

void FVMSolver::print_matrix(vector<vector<double>> *m){
    auto& mat = *m;
    int n = mat.size();

    cout << "[" << endl;
    for (int i = 0; i < n; i++)
    {
        cout << " [";
        for (int j = 0; j < n; j++)
        {
            cout << setw(8) << fixed << setprecision(3) << mat[i][j];
            if (j < n - 1)
                cout << ", ";
        }
        cout << "]" << (i < n - 1 ? "," : "") << endl;
    }
    cout << "]" << endl;
}

void FVMSolver::print_vector(vector<double> *v){
    auto& vet = *v;
    cout << " [";
    int n = vet.size();
    
    for(int i = 0; i < n; i++){
        cout << setw(8) << fixed << setprecision(3) << vet[i];
        if (i < n - 1)
                cout << ", ";
    }
    cout << "]" << endl;
}

void FVMSolver::assembly_A(){
    int ncells = mesh->get_num_cells();
    for(int i = 0; i < ncells; i++){
        int gcell = mesh->get_global_cell_id(i); 
        vector<int>& faceIds = mesh->get_cell(gcell)->get_face_ids();
        
        for(int j = 0; j < faceIds.size(); j++){
            int gface = faceIds[j];

            if(mesh->get_link_face_to_bface(gface) != -1){
                // Contribuição da face de contorno
                Node* centroid = mesh->get_centroid(i);
                Node* middleFace = mesh->get_middle_face(gface);

                double lb1 = middleFace->get_x() - centroid->get_x();
                double lb2 = middleFace->get_y() - centroid->get_y();

                tuple<double, double>& normal = mesh->get_normal(gface);

                double deltab = lb1 * get<0>(normal) + lb2 * get<1>(normal);

                A[i][i] += (gamma[i] * mesh->get_face_area(gface)) / deltab;
            }else{
                // Contribuição da face do interior
                pair<int,int>& idCellsShareFace = mesh->get_link_face_to_cell(gface);

                int ic1 = idCellsShareFace.first; 
                int ic2 = idCellsShareFace.second;
                
                int nb = ic1 == gcell ? ic2 : ic1; //pegando o vizinho do gcell
                int lnb = mesh->get_local_cell_id(nb);

                A[i][i] += (gamma[i] * mesh->get_face_area(gface)) / mesh->get_deltaf(gface); // diagonal
                A[i][lnb] = - (gamma[i] *  mesh->get_face_area(gface)) / this->mesh->get_deltaf(gface); // off-diagonal
            }
        }
    }
}

void FVMSolver::assembly_b(){
    int ncells = mesh->get_num_cells();
    for(int i = 0; i < ncells; i++){
        int gcell = mesh->get_global_cell_id(i);
        vector<int>& faceIds = mesh->get_cell(gcell)->get_face_ids();
        for(int j = 0; j < faceIds.size(); j++){
            int gface = faceIds[j];
            
            if(mesh->get_link_face_to_bface(gface) != -1) // face de contorno
                continue; // não faz nada

            /*Parte do skew... (REVER ISSO DEPOIS)*/
            pair<int, int>& idNodes = mesh->get_link_face_to_vertex(gface);
            int va = idNodes.first;
            int vb = idNodes.second;

            pair<int, int>& idCellsShareFace = mesh->get_link_face_to_cell(gface);
            int ic1 = idCellsShareFace.first; 
            int ic2 = idCellsShareFace.second;
            
            Node* centroid1 = this->mesh->get_centroid(i);
            Node* centroid2 = this->mesh->get_centroid(i);

            // calcula vetor Lb
            double Lb1 = centroid2->get_x() - centroid1->get_x();
            double Lb2 = centroid2->get_y() - centroid1->get_y();

            Node* vaNode = mesh->get_node(va);
            Node* vbNode = mesh->get_node(vb);

            // calcula vetor tangente a face (unitária)
            double tf1 = (vaNode->get_x() - vbNode->get_x()) / this->mesh->get_face_area(gface);
            double tf2 = (vaNode->get_y() - vbNode->get_y()) / this->mesh->get_face_area(gface);

            double dot_tf_Lb = tf1 * Lb1 + tf2 * Lb2;

            skew[i] += -(this->gamma[i] * dot_tf_Lb * (phiv(vaNode) - phiv(vbNode))) / mesh->get_deltaf(gface);
        }
    }
}

/*Computa a norma infinito da diferença entre dois vetores u1 e u2*/
double FVMSolver::max_norm_diff(vector<double>& u1, vector<double>& u2){
    double max_diff = 0.0;
    for (int i = 0; i < u1.size(); ++i) {
        max_diff = max(max_diff, fabs(u1[i] - u2[i]));
    }
    return max_diff;
}

void FVMSolver::save_solution(string filename){
    
    ofstream vtk_file(filename);

    if (!vtk_file.is_open()) {
        std::cerr << "ERROR: failed to create file." << filename << "'.\n";
        exit(1);
    }

    /*Informações padrão (REVER A QUESTÃO DA TEMPERATURA DEPOIS...)*/
    vtk_file << "# vtk DataFile Version 3.0\n";
    vtk_file << "Malha 2D com triângulos e temperatura\n";
    vtk_file << "ASCII\n";
    vtk_file << "DATASET UNSTRUCTURED_GRID\n\n";

    int nnodes = this->mesh->get_num_nodes();
    vtk_file << "POINTS " << nnodes << " double" << endl;
    vector<Node>& nodes = this->mesh->get_nodes();
    for(int i = 0; i < nnodes; i++){
        vtk_file << nodes[i].get_x() << " " << nodes[i].get_y() << " " << nodes[i].get_z() << endl;
    }
    vtk_file << endl;

    int ncells = this->mesh->get_num_cells();
    vtk_file << "CELLS " << ncells << " " <<  4 * ncells << endl;
    for(int i = 0; i < ncells; i++){
        int gcell = this->mesh->get_global_cell_id(i);
        Element* e = this->mesh->get_cell(gcell);
        vector<int>& nodeIds = e->get_nodes();
        vtk_file << "3" << " " << nodeIds[0] << " " << nodeIds[1] << " " << nodeIds[2] << endl;
    }
    vtk_file << endl;

    vtk_file << "CELL_TYPES " << ncells << endl;
    for(int i = 0; i < ncells; i++){
        vtk_file << "5" << endl;
    }
    vtk_file << endl;
    
    vtk_file << "CELL_DATA " << ncells << endl;
    vtk_file << "SCALARS Temperatura double 1\n";
    vtk_file << "LOOKUP_TABLE default\n";
    for(int i = 0; i < ncells; i++){
        vtk_file << this->u[i] << endl;
    }

    vtk_file.close();
}