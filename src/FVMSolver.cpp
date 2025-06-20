#include "../include/FVMSolver.h"

FVMSolver::FVMSolver(string filepath, BoundaryCondition *down, BoundaryCondition* right, BoundaryCondition* top, BoundaryCondition *left, double (*g)(double, double), double (*sourceTerm)(double, double)){
    this->mesh = new Mesh();
    mesh->read_mesh(filepath);

    /*Salva as condições usando a mesma orientação definida.*/
    this->boundaries.push_back(down);
    this->boundaries.push_back(right);
    this->boundaries.push_back(top);
    this->boundaries.push_back(left);

    mesh->pre_processing(&boundaries);

    this->gammafunc = g;
    this->sourcefunc = sourceTerm;

    int ncells = this->mesh->get_num_cells();
    int nfaces = mesh->get_num_faces();

    /*Inicialização das EDs associadas a resolução do problema...*/
    this->A = Eigen::MatrixXd(ncells, ncells); // ncells x ncells : 0's
    this->b = Eigen::VectorXd(ncells); // 1 x n : 0's
    this->u = Eigen::VectorXd(ncells); // 1 x n: 0's
    this->source = Eigen::VectorXd(ncells); // 1 x n: 0's
    this->gammaf = Eigen::VectorXd(nfaces); 
    this->skew = Eigen::VectorXd(ncells);

    // chama para fazer o pré-processamento
    this->pre_processing();
}

FVMSolver::~FVMSolver(){
}

void FVMSolver::pre_processing(){
    int ncells = mesh->get_num_cells();
    int nfaces = mesh->get_num_faces();

    int gcell;
    Element* cell;
    double source1, source2, source3;
    double gamma1, gamma2, gamma3;
    for(int i = 0; i < ncells; i++){ // itera pelos índices locais
        int gcell = mesh->get_global_cell_id(i);
        Element* cell = mesh->get_cell(gcell);
        vector<int>& faces = cell->get_face_ids();
        double sum = 0;
        for(int j = 0; j < faces.size(); j++){
            Node* midFace = mesh->get_middle_face(faces[j]);
            sum += this->sourcefunc(midFace->get_x(), midFace->get_y());
        }

        this->source[i] = sum / faces.size();
    }

    for(int i = 0; i < nfaces; i++){
        pair<int,int>& nodeIds = mesh->get_link_face_to_vertex(i);
        Node* n1 = mesh->get_node(nodeIds.first);
        Node* n2 = mesh->get_node(nodeIds.second);
        Node* middle = mesh->get_middle_face(i);

        double d1 = sqrt(pow(n1->get_x() - middle->get_x(), 2) + pow(n1->get_y() - middle->get_y(), 2));
        double d2 = sqrt(pow(n2->get_x() - middle->get_x(), 2) + pow(n2->get_y() - middle->get_y(), 2));

        double gamma1 = gammafunc(n1->get_x(), n1->get_y());
        double gamma2 = gammafunc(n2->get_x(), n2->get_y());

        this->gammaf[i] = (gamma1/d1 + gamma2/d2)/(1/d1 + 1/d2);
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
            cout << setw(8) << fixed << setprecision(8) << mat[i][j];
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
        cout << setw(8) << fixed << setprecision(8) << vet[i];
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

                A(i,i) += (gammaf[gface] * mesh->get_face_area(gface)) / deltab;
            }else{
                // Contribuição da face do interior
                pair<int,int>& idCellsShareFace = mesh->get_link_face_to_cell(gface);
                
                int ic1 = idCellsShareFace.first; 
                int ic2 = idCellsShareFace.second;
                
                int nb = ic1 == gcell ? ic2 : ic1; //pegando o vizinho do gcell
                int lnb = mesh->get_local_cell_id(nb);

                A(i,i) += (gammaf[gface] * mesh->get_face_area(gface)) / mesh->get_deltaf(gface); // diagonal
                A(i, lnb) = -(gammaf[gface] *  mesh->get_face_area(gface)) / this->mesh->get_deltaf(gface); // off-diagonal
            }
        }
    }
}

void FVMSolver::assembly_b(){
    int ncells = mesh->get_num_cells();
    for(int i = 0; i < ncells; i++){
        int gcell = mesh->get_global_cell_id(i);
        // vector<int>& faceIds = mesh->get_cell(gcell)->get_face_ids();
        // for(int j = 0; j < faceIds.size(); j++){
        //     int gface = faceIds[j];
            
        //     if(mesh->get_link_face_to_bface(gface) != -1) // face de contorno
        //         continue; // não faz nada

        //     /*Parte do skew... (REVER ISSO DEPOIS), pois por enquanto está tudo zerado e não atualiza.*/
        //     pair<int, int>& idNodes = mesh->get_link_face_to_vertex(gface);
        //     int va = idNodes.first;
        //     int vb = idNodes.second;

        //     pair<int, int>& idCellsShareFace = mesh->get_link_face_to_cell(gface);
        //     int ic1 = idCellsShareFace.first; 
        //     int ic2 = idCellsShareFace.second;
            
        //     Node* centroid1 = this->mesh->get_centroid(mesh->get_local_cell_id(ic1));
        //     Node* centroid2 = this->mesh->get_centroid(mesh->get_local_cell_id(ic2));

        //     // calcula vetor Lb
        //     double Lb1 = centroid2->get_x() - centroid1->get_x();
        //     double Lb2 = centroid2->get_y() - centroid1->get_y();

        //     Node* vaNode = mesh->get_node(va);
        //     Node* vbNode = mesh->get_node(vb);

        //     // calcula vetor tangente a face (unitária)
        //     double tf1 = (vaNode->get_x() - vbNode->get_x()) / this->mesh->get_face_area(gface);
        //     double tf2 = (vaNode->get_y() - vbNode->get_y()) / this->mesh->get_face_area(gface);

        //     double dot_tf_Lb = tf1 * Lb1 + tf2 * Lb2;

        //     skew[i] += -(this->gammaf[gface] * dot_tf_Lb * (phiv(vaNode) - phiv(vbNode))) / mesh->get_deltaf(gface);
        // }
        b[i] = - source[i] * mesh->get_volume(i);
    }
    this->apply_boundaries();
}

void FVMSolver::apply_boundaries(){
    int nbfaces = this->mesh->get_num_boundary_faces();
    for(int i = 0; i < nbfaces; i++){ //PARA TODA FACE DE CONTORNO NUMERADA LOCALMENTE
        int gface = this->mesh->get_link_bface_to_face(i); //PEGA A NUMERACAO GLOBAL
        pair<int,int>& cell = this->mesh->get_link_face_to_cell(gface); //PEGA AS CELULAS QUE FAZEM INTERSEÇÃO NAQUELA FACE
        int gcell = cell.first != -1 ? cell.first : cell.second;
        int lcell = mesh->get_local_cell_id(gcell);

        Node* centroid = mesh->get_centroid(lcell);
        Node* middleFace = mesh->get_middle_face(gface);

        double lb1 = middleFace->get_x() - centroid->get_x();
        double lb2 = middleFace->get_y() - centroid->get_y();

        tuple<double, double>& normal = mesh->get_normal(gface);

        double deltab = lb1 * get<0>(normal) + lb2 * get<1>(normal);

        b[lcell] += gammaf[gface] * mesh->get_face_area(gface) * mesh->get_phib(i) / deltab;
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

void FVMSolver::iterative_solver(){
    Eigen::LDLT<Eigen::MatrixXd> ldlt(A);
    if(ldlt.info() != Eigen::Success) {
        std::cerr << "Erro na decomposição LDLT!" << std::endl;
        exit(-1);
    }
    this->u = ldlt.solve(b);

    if(ldlt.info() != Eigen::Success) {
        std::cerr << "Erro na resolução do sistema!" << std::endl;
        exit(-1);
    }

    // cout << "Simetrica: " <<  A.isApprox(A.transpose()) << endl;
    // u = A.LDLT().solve(b);
    // std::cout.precision(8);
    // std::cout << std::fixed;

    // Imprima a solução x
    // cout << "Determinante: " << A.determinant() << endl;
    // std::cout << "A solucao u é:" << std::endl;
    std::cout << u << std::endl;
}

double FVMSolver::compute_error(double (*exact)(double, double)){
    Eigen::VectorXd exact_vect(mesh->get_num_cells());
    for(int i = 0; i < mesh->get_num_cells(); i++){ //para cada célula
        int gcell = mesh->get_global_cell_id(i);
        Element* cell = mesh->get_cell(gcell);
        vector<int>& faces = cell->get_face_ids();
        double sum = 0;
        for(int j = 0; j < faces.size(); j++){
            Node* midFace = mesh->get_middle_face(faces[j]);
            sum += exact(midFace->get_x(), midFace->get_y());
        }
        exact_vect[i] = sum / faces.size();
    }

    /*calcula norma do máximo*/
    double max_diff = 0.0;
    for (int i = 0; i < mesh->get_num_cells(); ++i) {
        max_diff = max(max_diff, fabs(u[i] - exact_vect[i]));
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
    vtk_file << "CELLS " << ncells << " ";
    int sum = 0;
    for(int i = 0; i < ncells; i++)
    {
        int gcell = this->mesh->get_global_cell_id(i);
        Element* e = this->mesh->get_cell(gcell);
        if(e->get_element_type() == 2)
            sum += 4;
        else if(e->get_element_type() == 3)
            sum += 5;
    }
    vtk_file << sum << endl;
    for(int i = 0; i < ncells; i++){
        int gcell = this->mesh->get_global_cell_id(i);
        Element* e = this->mesh->get_cell(gcell);
        vector<int>& nodeIds = e->get_nodes();
        vtk_file << nodeIds.size(); 
        for(int j = 0; j < nodeIds.size(); j++)
            vtk_file << " " << nodeIds[j];
        vtk_file << endl;
    }
    vtk_file << endl;

    vtk_file << "CELL_TYPES " << ncells << endl;
    for(int i = 0; i < ncells; i++){
        int gcell = mesh->get_global_cell_id(i);
        Element* cell = mesh->get_cell(gcell);
        if(cell->get_element_type() == 2)
            vtk_file << "5" << endl;
        else if(cell->get_element_type() == 3)
            vtk_file << "9" << endl;
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

bool FVMSolver::is_simetric(){
    for(int i = 0, j = this->u.size() - 1; i < j; i++, j--)
    {
        if(u[i] - u[j] > 1e-2)
            return false;
    }
    return true;
}