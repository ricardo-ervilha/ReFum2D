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
    this->A = arma::mat(ncells, ncells); // ncells x ncells : 0's
    this->b = arma::vec(ncells); // 1 x n : 0's
    this->u = arma::vec(ncells); // 1 x n: 0's
    this->source = arma::vec(ncells); // 1 x n: 0's
    this->gammaf = arma::vec(nfaces); 
    this->b_with_cd = arma::vec(ncells);
    
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
        vector<int>& normalSigns = mesh->get_cell(gcell)->get_normal_sign();
        for(int j = 0; j < faceIds.size(); j++){
            int gface = faceIds[j];

            if(mesh->get_link_face_to_bface(gface) != -1){
                // Contribuição da face de contorno
                Node* centroid = mesh->get_centroid(i);
                Node* middleFace = mesh->get_middle_face(gface);
                
                /*cálculo da n1*/
                tuple<double, double>& normal = mesh->get_normal(gface);
                
                //correção do sinal da normal para apontar para fora de P corretamente.
                double nx = get<0>(normal)*normalSigns[j];
                double ny = get<1>(normal)*normalSigns[j];
                
                double dpnx = middleFace->get_x() - centroid->get_x(); 
                double dpny = middleFace->get_y() - centroid->get_y();

                double normnfsquare = nx*nx + ny*ny;
                
                double nfdotdpn = nx*dpnx + ny*dpny;
                
                double n1x = dpnx * (normnfsquare/nfdotdpn);
                double n1y = dpny * (normnfsquare/nfdotdpn);

                double normn1 = n1x*n1x + n1y*n1y;
                /*fim cálculo da n1*/

                // distância entre o centroide de P e a face de contorno.
                double deltab = sqrt(pow(centroid->get_x() - middleFace->get_x(), 2) + pow(centroid->get_y() - middleFace->get_y(), 2));

                A(i,i) -= (gammaf[gface] * mesh->get_face_area(gface) * normn1) / deltab;
            }else{
                // Contribuição de uma face que NÃO é face de contorno.
                pair<int,int>& idCellsShareFace = mesh->get_link_face_to_cell(gface);
                
                int ic1 = idCellsShareFace.first; 
                int ic2 = idCellsShareFace.second;
                
                int nb = ic1 == gcell ? ic2 : ic1; //pegando o vizinho do gcell
                int lnb = mesh->get_local_cell_id(nb);
                
                /*cálculo de |n1|*/
                Node *centroidP = mesh->get_centroid(i);
                Node *centroidN = mesh->get_centroid(lnb);
                tuple<double, double>& normal = mesh->get_normal(gface);
                
                //correção do sinal da normal para apontar para fora de P corretamente.
                double nx = get<0>(normal)*normalSigns[j];
                double ny = get<1>(normal)*normalSigns[j];

                double dpnx = centroidN->get_x() - centroidP->get_x();
                double dpny = centroidN->get_y() - centroidP->get_y();
                
                double normnfsquare = nx*nx + ny*ny;
                double nfdotdpn = nx*dpnx + ny*dpny;
                
                double n1x = dpnx * (normnfsquare/nfdotdpn);
                double n1y = dpny * (normnfsquare/nfdotdpn);

                double normn1 = n1x*n1x + n1y*n1y;
                /*fim do cálculo de |n1|*/

                A(i,i) -= (gammaf[gface] * mesh->get_face_area(gface) * normn1) / mesh->get_deltaf(gface); // diagonal
                A(i, lnb) = (gammaf[gface] *  mesh->get_face_area(gface) * normn1) / this->mesh->get_deltaf(gface); // off-diagonal
            }
        }
    }
}

void FVMSolver::assembly_b(){
    int ncells = mesh->get_num_cells();
    for(int i = 0; i < ncells; i++){
        int gcell = mesh->get_global_cell_id(i);
        
        b[i] = source[i] * mesh->get_volume(i);
    }
    this->apply_boundaries();
}

/*Uso da reconstrução com least squares, e interpolação da phif usando os valores dos centroides.*/
/*Least squares irá resolver Mx=y*/
void FVMSolver::compute_gradients(){
    int ncells = mesh->get_num_cells();
    this->gradients = vector<std::pair<double, double>>(ncells, std::make_pair(0.0, 0.0));

    //construirá o vetor b_with_cd para ter a difusão cruzada.
    for(int i = 0; i < ncells; i++){ // para cada célula
        int gcell = mesh->get_global_cell_id(i);
        Element* cell = mesh->get_cell(gcell);
        vector<int> faceIds = cell->get_face_ids();
        
        arma::mat M(faceIds.size(), 2);
        arma::vec y(faceIds.size());
        
        for(int j = 0; j < faceIds.size(); j++){ // para cada face
            int gface = faceIds[j];
            if(mesh->get_link_face_to_bface(gface) != -1){
                // Contribuição da face de contorno
                Node* centroid = mesh->get_centroid(i);
                Node* middleFace = mesh->get_middle_face(gface);
                
                double dx = middleFace->get_x() - centroid->get_x();
                double dy = middleFace->get_y() - centroid->get_y(); 

                M(j,0) = dx;
                M(j, 1) = dy;
                y[j] = -u[i]; // u_B (0 por enquanto !) - u_P  
            }else{
                pair<int,int>& idCellsShareFace = mesh->get_link_face_to_cell(gface);
                    
                int ic1 = idCellsShareFace.first; 
                int ic2 = idCellsShareFace.second;
                
                int nb = ic1 == gcell ? ic2 : ic1; //pegando o vizinho do gcell
                int lnb = mesh->get_local_cell_id(nb);
                
                /*cálculo de |n1|*/
                Node *centroidP = mesh->get_centroid(i);
                Node *centroidN = mesh->get_centroid(lnb);

                double dx = centroidN->get_x() - centroidP->get_x();
                double dy = centroidN->get_y() - centroidP->get_y();

                M(j,0) = dx;
                M(j, 1) = dy;
                y[j] = u[lnb] - u[i]; //u_N - u_P  
            }
        }

        // M x = y
        arma::vec x = arma::solve(M, y);
        this->gradients[i] = make_pair(x[0], x[1]);
    }
}

void FVMSolver::compute_cross_diffusion(){
    int ncells = mesh->get_num_cells();
    b_with_cd.zeros();
    for(int i = 0; i < ncells; i++){
        int gcell = mesh->get_global_cell_id(i);
        Element* cell = mesh->get_cell(gcell);
        vector<int> faceIds = cell->get_face_ids();
        vector<int>& normalSigns = mesh->get_cell(gcell)->get_normal_sign();
        for(int j = 0; j < faceIds.size(); j++){ // para cada face
            int gface = faceIds[j];
            if(mesh->get_link_face_to_bface(gface) != -1){
                //face de contorno
                Node* centroid = mesh->get_centroid(i);
                Node* middleFace = mesh->get_middle_face(gface);
                
                /*cálculo da n1*/
                tuple<double, double>& normal = mesh->get_normal(gface);
                
                //correção do sinal da normal para apontar para fora de P corretamente.
                double nx = get<0>(normal)*normalSigns[j];
                double ny = get<1>(normal)*normalSigns[j];
                
                double dpnx = middleFace->get_x() - centroid->get_x(); 
                double dpny = middleFace->get_y() - centroid->get_y();

                double normnfsquare = nx*nx + ny*ny;
                
                double nfdotdpn = nx*dpnx + ny*dpny;
                
                double n1x = dpnx * (normnfsquare/nfdotdpn);
                double n1y = dpny * (normnfsquare/nfdotdpn);

                double n2x = nx - n1x;
                double n2y = ny - n1y;

                double graddotnormal = gradients[i].first * n2x + gradients[i].second * n2y;

                this->b_with_cd[i] += (gammaf[gface] * mesh->get_face_area(gface) * graddotnormal);
            }else{
                pair<int,int>& idCellsShareFace = mesh->get_link_face_to_cell(gface);
                Node* middleFace = mesh->get_middle_face(gface);
                
                int ic1 = idCellsShareFace.first; 
                int ic2 = idCellsShareFace.second;
                
                int nb = ic1 == gcell ? ic2 : ic1; //pegando o vizinho do gcell
                int lnb = mesh->get_local_cell_id(nb);
                
                /*cálculo de |n1|*/
                Node *centroidP = mesh->get_centroid(i);
                Node *centroidN = mesh->get_centroid(lnb);
                tuple<double, double>& normal = mesh->get_normal(gface);
                
                //correção do sinal da normal para apontar para fora de P corretamente.
                double nx = get<0>(normal)*normalSigns[j];
                double ny = get<1>(normal)*normalSigns[j];

                double dpnx = centroidN->get_x() - centroidP->get_x();
                double dpny = centroidN->get_y() - centroidP->get_y();
                
                double normnfsquare = nx*nx + ny*ny;
                double nfdotdpn = nx*dpnx +ny*dpny;
                
                double n1x = dpnx * (normnfsquare/nfdotdpn);
                double n1y = dpny * (normnfsquare/nfdotdpn);

                double n2x = nx - n1x;
                double n2y = ny - n1y;

                double d1 = sqrt(pow(middleFace->get_x() - centroidN->get_x(),2) + pow(middleFace->get_y() - centroidN->get_y(),2));
                double d2 = sqrt(pow(middleFace->get_x() - centroidP->get_x(),2) + pow(middleFace->get_y() - centroidP->get_y(),2));

                double gradx = (gradients[i].first*d1 + gradients[lnb].first*d2)/(d1+d2);
                double grady = (gradients[i].second*d1 + gradients[lnb].second*d2)/(d1+d2);

                double graddotnormal = gradx * n2x + grady * n2y;

                this->b_with_cd[i] += (gammaf[gface] * mesh->get_face_area(gface) * graddotnormal);
            }
        }
        this->b_with_cd[i] = this->b[i] - this->b_with_cd[i];
    }
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

    //CORRIGIR ESSA FUNÇÃO POIS ELA ESTÁ DESATUALIZADA
}

void FVMSolver::solve_system(){
    for(int k = 0; k < 100; k++){
        
        this->compute_gradients();
        this->compute_cross_diffusion();
        this->u = arma::solve(A,b);
    }

    cout << "\nSolução obtida:\n";
    cout << u << endl;
}

double FVMSolver::compute_error(double (*exact)(double, double)){
    arma::vec exact_vect(mesh->get_num_cells());
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
    double norml2 = 0;
    double sumAreas = 0;
    for (int i = 0; i < mesh->get_num_cells(); ++i) {
        max_diff = max(max_diff, fabs(u[i] - exact_vect[i]));
        norml2 += pow(u[i] - exact_vect[i],2) * mesh->get_volume(i);
        sumAreas +=  mesh->get_volume(i);
    }
    norml2 = sqrt(norml2/sumAreas);
    cout << "\nNorma L2: " << norml2 << endl;
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