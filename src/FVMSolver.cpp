#include "../include/FVMSolver.h"

FVMSolver::FVMSolver(Mesh *mesh, BoundaryCondition *down, BoundaryCondition* right, BoundaryCondition* top, BoundaryCondition *left, double (*g)(double, double), double(*rho)(double,double), pair<double,double>(*U)(double,double), double (*sourceTerm)(double, double)){
    this->mesh = mesh;

    /*Salva as condições usando a mesma orientação definida.*/
    this->boundaries.push_back(down);
    this->boundaries.push_back(right);
    this->boundaries.push_back(top);
    this->boundaries.push_back(left);

    this->gammafunc = g;
    this->sourcefunc = sourceTerm;
    this->rhofunc = rho;
    this->ufunc = U;

    int ncells = this->mesh->get_ncells();
    int nfaces = mesh->get_nedges();

    /*Inicialização das EDs associadas a resolução do problema...*/
    this->A = arma::mat(ncells, ncells); 
    this->b = arma::vec(ncells);
    this->u = arma::vec(ncells); 
    this->source = arma::vec(ncells); 
    this->gammaf = arma::vec(nfaces); 
    this->rhof = arma::vec(nfaces); 
    this->uf = arma::mat(nfaces,2); 
    this->b_with_cd = arma::vec(ncells);
    
    // chama para fazer o pré-processamento
    this->pre_processing();
    // chama também para encontrar os phibs.
    this->apply_boundaries_in_edges_from_mesh();
}

FVMSolver::~FVMSolver(){
}

void FVMSolver::pre_processing(){
    int ncells = mesh->get_ncells();
    int nfaces = mesh->get_nedges();

    vector<Cell*>& cells = mesh->get_cells();
    vector<Edge*>& edges = mesh->get_edges();

    for(int i = 0; i < ncells; i++){ // itera pelos índices locais
        pair<double,double>& centroid = cells[i]->get_centroid();
        this->source[i] = this->sourcefunc(centroid.first, centroid.second);
    }
    
    for(int i = 0; i < nfaces; i++){
        if(edges[i]->is_boundary_face()){
            //face de contorno
            pair<double,double>& midface = edges[i]->get_middle();
            gammaf[i] = gammafunc(midface.first, midface.second);
            rhof[i] = gammafunc(midface.first, midface.second);
            pair<double,double> u = ufunc(midface.first, midface.second);
            uf(i,0) = u.first; uf(i,1) = u.second;
        }else{
            //não é face de contorno
            pair<double,double>& midface = edges[i]->get_middle();
            pair<int,int>& idCellsShareFace = edges[i]->get_link_face_to_cell();
              
            pair<double,double>& centroidP = cells[idCellsShareFace.first]->get_centroid();
            pair<double,double>& centroidN = cells[idCellsShareFace.second]->get_centroid();

            double d1 = distance(midface.first, midface.second, centroidP.first, centroidP.second);
            double d2 = distance(midface.first, midface.second, centroidN.first, centroidN.second);
            double totaldist = d1+d2;

            gammaf[i] = (gammafunc(centroidP.first, centroidP.second)*d2 + gammafunc(centroidN.first, centroidN.second)*d1)/totaldist;
            rhof[i] = (rhofunc(centroidP.first, centroidP.second)*d2 + rhofunc(centroidN.first, centroidN.second)*d1)/totaldist;
            uf(i,0) = (ufunc(centroidP.first, centroidP.second).first*d2 + ufunc(centroidN.first, centroidN.second).first*d1)/totaldist;
            uf(i,1) = (ufunc(centroidP.first, centroidP.second).second*d2 + ufunc(centroidN.first, centroidN.second).second*d1)/totaldist;
        }
    }
}

void FVMSolver::apply_boundaries_in_edges_from_mesh(){
    vector<Edge*>& edges = mesh->get_edges();
    for(int i = 0; i < edges.size(); i++){
        if(edges[i]->is_boundary_face()){
            // é face de contorno
            pair<double,double> midFace = edges[i]->get_middle();
            if(edges[i]->from->y == mesh->get_ymin() && edges[i]->to->y == mesh->get_ymin()){
                // down boundary
                edges[i]->set_phib(this->boundaries[0]->apply(midFace.first, midFace.second));
            } else if(edges[i]->from->x == mesh->get_xmax() && edges[i]->to->x == mesh->get_xmax()){
                // right boundary
                edges[i]->set_phib(this->boundaries[1]->apply(midFace.first, midFace.second));
            } else if(edges[i]->from->y == mesh->get_ymax() && edges[i]->to->y == mesh->get_ymax()){
                // top boundary
                edges[i]->set_phib(this->boundaries[2]->apply(midFace.first, midFace.second));
            } else if(edges[i]->from->x == mesh->get_xmin() && edges[i]->to->x == mesh->get_xmin()){
                //left boundary
                edges[i]->set_phib(this->boundaries[3]->apply(midFace.first, midFace.second));
            }
        }
    }
}

pair<double,double> compute_n1(pair<double,double>& p, pair<double,double>& n, pair<double,double>& normal){
    pair<double,double> dpn = make_pair(n.first - p.first, n.second - p.second);
    double norm_normal_square = normal.first*normal.first + normal.second*normal.second; // |nf|^2
    double nfdotdpn = normal.first*dpn.first + normal.second*dpn.second; // nf . dpn
    
    double n1x = dpn.first * (norm_normal_square/nfdotdpn);
    double n1y = dpn.second * (norm_normal_square/nfdotdpn);

    return make_pair(n1x, n1y);
}

void FVMSolver::assembly_A(){
    vector<Cell*>& cells = mesh->get_cells();
    for(int i = 0; i < cells.size(); i++){
        vector<Edge*> facesOfCell = cells[i]->get_edges();
        vector<int>& nsigns = cells[i]->get_nsigns();
        
        for(int j = 0; j < facesOfCell.size(); j++){
            int gface = facesOfCell[j]->id;
            if(facesOfCell[j]->is_boundary_face()){
                // Contribuição da face de contorno
                pair<double,double>& centroid = cells[i]->get_centroid();
                pair<double,double>& middleFace = facesOfCell[j]->get_middle();
                pair<double, double>& normal = facesOfCell[j]->get_normal();
                
                // correção da normal
                pair<double,double> normal_corrected = make_pair(normal.first *  nsigns[j], normal.second * nsigns[j]);
                
                pair<double,double> n1 = compute_n1(centroid, middleFace, normal_corrected);
                double normn1 = sqrt(n1.first * n1.first + n1.second * n1.second); // |n1|

                // distância entre o centroide de P e a face de contorno.
                double deltab = distance(centroid.first, centroid.second, middleFace.first, middleFace.second);
                
                /*Difusão.*/
                // - k_B * A_B * |n1| / deltab
                A(i,i) += -(gammaf[gface] * facesOfCell[j]->get_length() * normn1) / deltab;
            }else{
                pair<int,int>& idCellsShareFace = facesOfCell[j]->get_link_face_to_cell();
                
                int ic1 = idCellsShareFace.first; 
                int ic2 = idCellsShareFace.second;
                
                int nb = ic1 == i ? ic2 : ic1; //pegando o vizinho da célula P
                
                pair<double,double>& centroidP = cells[i]->get_centroid();
                pair<double,double>& centroidN = cells[nb]->get_centroid();
                pair<double, double>& normal = facesOfCell[j]->get_normal();
                pair<double,double> normal_corrected = make_pair(normal.first *  nsigns[j], normal.second * nsigns[j]);
                
                pair<double,double> n1 = compute_n1(centroidP, centroidN, normal_corrected);
                double normn1 = sqrt(n1.first * n1.first + n1.second * n1.second); // |n1|
                
                /*Calculo da distancia de P até center(gface) & N até center(gface)*/
                pair<double,double>& midFace = facesOfCell[j]->get_middle();
                double d1 = distance(centroidP.first, centroidP.second, midFace.first, midFace.second);
                double d2 = distance(centroidN.first, centroidN.second, midFace.first, midFace.second);
                
                /*Difusão*/
                // na diagonal: -k_f * A_f * |n1| / delta_f
                double coeff = (gammaf[gface] * facesOfCell[j]->get_length() * normn1) / facesOfCell[j]->get_df();
                A(i,i) += -coeff;

                // fora da diagonal: +k_f * A_f * |n1| / delta_f
                A(i, nb) += coeff; 
            }
        }
    }
}

void FVMSolver::assembly_b(){
    
    vector<Cell*>& cells = mesh->get_cells();
    for(int i = 0; i < cells.size(); i++){
        vector<Edge*> facesOfCell = cells[i]->get_edges();
        vector<int>& nsigns = cells[i]->get_nsigns();
        for(int j = 0; j < facesOfCell.size(); j++){
            int gface = facesOfCell[j]->id;
            if(facesOfCell[j]->is_boundary_face()){
                // caso seja face de contorno
                pair<double,double>& centroid = cells[i]->get_centroid();
                pair<double,double>& middleFace = facesOfCell[j]->get_middle();
                pair<double, double>& normal = facesOfCell[j]->get_normal();
                
                // correção da normal
                pair<double,double> normal_corrected = make_pair(normal.first *  nsigns[j], normal.second * nsigns[j]);
                
                pair<double,double> n1 = compute_n1(centroid, middleFace, normal_corrected);
                double normn1 = sqrt(n1.first * n1.first + n1.second * n1.second); // |n1|

                // distância entre o centroide de P e a face de contorno.
                double deltab = distance(centroid.first, centroid.second, middleFace.first, middleFace.second);

                // no b da celula i, será: - k_B * A_B * |n1| * phi_B / delta_B
                b[i] = - (gammaf[gface] * facesOfCell[j]->get_length() * normn1 * facesOfCell[j]->get_phib()) / deltab;
            }
        }
        // acrescento agora o termo fonte: S_P * V_P.
        b[i] += this->source[i] * cells[i]->get_area();
    }

    /*
        b possuirá a contribuição fixa.
        falta ainda a difusão cruzada, que será computada usando um vetor auxiliar.
    */
}

// /*Uso da reconstrução com least squares, e interpolação da phif usando os valores dos centroides.*/
// /*Least squares irá resolver Mx=y*/
// void FVMSolver::compute_gradients(){
//     int ncells = mesh->get_num_cells();
//     this->gradients = vector<std::pair<double, double>>(ncells, std::make_pair(0.0, 0.0));

//     //construirá o vetor b_with_cd para ter a difusão cruzada.
//     for(int i = 0; i < ncells; i++){ // para cada célula
//         int gcell = mesh->get_global_cell_id(i);
//         Element* cell = mesh->get_cell(gcell);
//         vector<int> faceIds = cell->get_face_ids();
        
//         arma::mat M(faceIds.size(), 2);
//         arma::vec y(faceIds.size());
        
//         for(int j = 0; j < faceIds.size(); j++){ // para cada face
//             int gface = faceIds[j];
//             if(mesh->get_link_face_to_bface(gface) != -1){
//                 // Contribuição da face de contorno
//                 Node* centroid = mesh->get_centroid(i);
//                 Node* middleFace = mesh->get_middle_face(gface);
                
//                 double dx = middleFace->get_x() - centroid->get_x();
//                 double dy = middleFace->get_y() - centroid->get_y(); 

//                 M(j,0) = dx;
//                 M(j, 1) = dy;
//                 y[j] = -u[i]; // u_B (0 por enquanto !) - u_P  
//             }else{
//                 pair<int,int>& idCellsShareFace = mesh->get_link_face_to_cell(gface);
                    
//                 int ic1 = idCellsShareFace.first; 
//                 int ic2 = idCellsShareFace.second;
                
//                 int nb = ic1 == gcell ? ic2 : ic1; //pegando o vizinho do gcell
//                 int lnb = mesh->get_local_cell_id(nb);
                
//                 /*cálculo de |n1|*/
//                 Node *centroidP = mesh->get_centroid(i);
//                 Node *centroidN = mesh->get_centroid(lnb);

//                 double dx = centroidN->get_x() - centroidP->get_x();
//                 double dy = centroidN->get_y() - centroidP->get_y();

//                 M(j,0) = dx;
//                 M(j, 1) = dy;
//                 y[j] = u[lnb] - u[i]; //u_N - u_P  
//             }
//         }

//         // M x = y
//         arma::vec x = arma::solve(M, y);
//         this->gradients[i] = make_pair(x[0], x[1]);
//     }
// }

// void FVMSolver::compute_cross_diffusion(){
//     int ncells = mesh->get_num_cells();
//     b_with_cd.zeros();
//     for(int i = 0; i < ncells; i++){
//         int gcell = mesh->get_global_cell_id(i);
//         Element* cell = mesh->get_cell(gcell);
//         vector<int> faceIds = cell->get_face_ids();
//         vector<int>& normalSigns = mesh->get_cell(gcell)->get_normal_sign();
//         for(int j = 0; j < faceIds.size(); j++){ // para cada face
//             int gface = faceIds[j];
//             if(mesh->get_link_face_to_bface(gface) != -1){
//                 //face de contorno
//                 Node* centroid = mesh->get_centroid(i);
//                 Node* middleFace = mesh->get_middle_face(gface);
                
//                 /*cálculo da n1*/
//                 tuple<double, double>& normal = mesh->get_normal(gface);
                
//                 //correção do sinal da normal para apontar para fora de P corretamente.
//                 double nx = get<0>(normal)*normalSigns[j];
//                 double ny = get<1>(normal)*normalSigns[j];
                
//                 double dpnx = middleFace->get_x() - centroid->get_x(); 
//                 double dpny = middleFace->get_y() - centroid->get_y();

//                 double normnfsquare = nx*nx + ny*ny;
                
//                 double nfdotdpn = nx*dpnx + ny*dpny;
                
//                 double n1x = dpnx * (normnfsquare/nfdotdpn);
//                 double n1y = dpny * (normnfsquare/nfdotdpn);

//                 double n2x = nx - n1x;
//                 double n2y = ny - n1y;

//                 double graddotnormal = gradients[i].first * n2x + gradients[i].second * n2y;

//                 this->b_with_cd[i] += (gammaf[gface] * mesh->get_face_area(gface) * graddotnormal);
//             }else{
//                 pair<int,int>& idCellsShareFace = mesh->get_link_face_to_cell(gface);
//                 Node* middleFace = mesh->get_middle_face(gface);
                
//                 int ic1 = idCellsShareFace.first; 
//                 int ic2 = idCellsShareFace.second;
                
//                 int nb = ic1 == gcell ? ic2 : ic1; //pegando o vizinho do gcell
//                 int lnb = mesh->get_local_cell_id(nb);
                
//                 /*cálculo de |n1|*/
//                 Node *centroidP = mesh->get_centroid(i);
//                 Node *centroidN = mesh->get_centroid(lnb);
//                 tuple<double, double>& normal = mesh->get_normal(gface);
                
//                 //correção do sinal da normal para apontar para fora de P corretamente.
//                 double nx = get<0>(normal)*normalSigns[j];
//                 double ny = get<1>(normal)*normalSigns[j];

//                 double dpnx = centroidN->get_x() - centroidP->get_x();
//                 double dpny = centroidN->get_y() - centroidP->get_y();
                
//                 double normnfsquare = nx*nx + ny*ny;
//                 double nfdotdpn = nx*dpnx +ny*dpny;
                
//                 double n1x = dpnx * (normnfsquare/nfdotdpn);
//                 double n1y = dpny * (normnfsquare/nfdotdpn);

//                 double n2x = nx - n1x;
//                 double n2y = ny - n1y;

//                 double d1 = sqrt(pow(middleFace->get_x() - centroidN->get_x(),2) + pow(middleFace->get_y() - centroidN->get_y(),2));
//                 double d2 = sqrt(pow(middleFace->get_x() - centroidP->get_x(),2) + pow(middleFace->get_y() - centroidP->get_y(),2));

//                 double gradx = (gradients[i].first*d1 + gradients[lnb].first*d2)/(d1+d2);
//                 double grady = (gradients[i].second*d1 + gradients[lnb].second*d2)/(d1+d2);

//                 double graddotnormal = gradx * n2x + grady * n2y;

//                 this->b_with_cd[i] += (gammaf[gface] * mesh->get_face_area(gface) * graddotnormal);
//             }
//         }
//         this->b_with_cd[i] = this->b[i] - this->b_with_cd[i];
//     }
// }

void FVMSolver::solve_system(){
   
    this->u = arma::solve(A,b);

    cout << "\nSolução obtida:\n";
    cout << u << endl;
}

void FVMSolver::compute_error(double (*exact)(double, double)){
    arma::vec exact_vect(mesh->get_ncells());
    vector<Cell*> cells = mesh->get_cells();
    for(int i = 0; i < cells.size(); i++){ //para cada célula
        pair<double,double>& centroid = cells[i]->get_centroid(); // obtém o centroide da celula
        exact_vect[i] = exact(centroid.first, centroid.second);
    }

    double norml2 = 0;
    double sumAreas = 0;
    for (int i = 0; i < cells.size(); ++i) {
        norml2 += (u[i] - exact_vect[i]) * (u[i] - exact_vect[i]) * cells[i]->get_area();
        sumAreas +=  cells[i]->get_area();
    }
    norml2 = sqrt(norml2/sumAreas);
    cout << "\nNorma L2: " << norml2 << endl;
}

// void FVMSolver::save_solution(string filename){
    
//     ofstream vtk_file(filename);

//     if (!vtk_file.is_open()) {
//         std::cerr << "ERROR: failed to create file." << filename << "'.\n";
//         exit(1);
//     }

//     /*Informações padrão (REVER A QUESTÃO DA TEMPERATURA DEPOIS...)*/
//     vtk_file << "# vtk DataFile Version 3.0\n";
//     vtk_file << "Malha 2D com triângulos e temperatura\n";
//     vtk_file << "ASCII\n";
//     vtk_file << "DATASET UNSTRUCTURED_GRID\n\n";

//     int nnodes = this->mesh->get_num_nodes();
//     vtk_file << "POINTS " << nnodes << " double" << endl;
//     vector<Node>& nodes = this->mesh->get_nodes();
//     for(int i = 0; i < nnodes; i++){
//         vtk_file << nodes[i].get_x() << " " << nodes[i].get_y() << " " << nodes[i].get_z() << endl;
//     }
//     vtk_file << endl;

//     int ncells = this->mesh->get_num_cells();
//     vtk_file << "CELLS " << ncells << " ";
//     int sum = 0;
//     for(int i = 0; i < ncells; i++)
//     {
//         int gcell = this->mesh->get_global_cell_id(i);
//         Element* e = this->mesh->get_cell(gcell);
//         if(e->get_element_type() == 2)
//             sum += 4;
//         else if(e->get_element_type() == 3)
//             sum += 5;
//     }
//     vtk_file << sum << endl;
//     for(int i = 0; i < ncells; i++){
//         int gcell = this->mesh->get_global_cell_id(i);
//         Element* e = this->mesh->get_cell(gcell);
//         vector<int>& nodeIds = e->get_nodes();
//         vtk_file << nodeIds.size(); 
//         for(int j = 0; j < nodeIds.size(); j++)
//             vtk_file << " " << nodeIds[j];
//         vtk_file << endl;
//     }
//     vtk_file << endl;

//     vtk_file << "CELL_TYPES " << ncells << endl;
//     for(int i = 0; i < ncells; i++){
//         int gcell = mesh->get_global_cell_id(i);
//         Element* cell = mesh->get_cell(gcell);
//         if(cell->get_element_type() == 2)
//             vtk_file << "5" << endl;
//         else if(cell->get_element_type() == 3)
//             vtk_file << "9" << endl;
//     }
//     vtk_file << endl;
    
//     vtk_file << "CELL_DATA " << ncells << endl;
//     vtk_file << "SCALARS Temperatura double 1\n";
//     vtk_file << "LOOKUP_TABLE default\n";
//     for(int i = 0; i < ncells; i++){
//         vtk_file << this->u[i] << endl;
//     }

//     vtk_file.close();
// }