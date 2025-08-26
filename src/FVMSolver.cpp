#include "../include/FVMSolver.h"

FVMSolver::FVMSolver(Mesh *mesh, BoundaryCondition *bc1, BoundaryCondition* bc2, BoundaryCondition* bc3, BoundaryCondition *bc4, double (*g)(double, double), double(*rho)(double,double), pair<double,double>(*U)(double,double), double (*sourceTerm)(double, double)){
    this->mesh = mesh;

    /*Salva as condições usando a mesma orientação definida.*/
    this->boundaries.emplace(bc1->get_location(), bc1);
    this->boundaries.emplace(bc2->get_location(), bc2);
    this->boundaries.emplace(bc3->get_location(), bc3);
    this->boundaries.emplace(bc4->get_location(), bc4);

    BoundaryCondition::apply_mins_maxs(mesh->get_xmin(), mesh->get_xmax(), mesh->get_ymin(), mesh->get_ymax());

    this->gammafunc = g;
    this->sourcefunc = sourceTerm;
    this->rhofunc = rho;
    this->ufunc = U;

    int ncells = this->mesh->get_ncells();
    int nfaces = mesh->get_nedges();

    /*Inicialização das ED's associadas a resolução do problema...*/
    // ! Matriz que guarda os coeficientes do sistema
    this->A = arma::sp_mat(ncells, ncells); // * ncells x ncells
    // ! Vetor coluna b, do sistema A u = b
    this->b = arma::vec(ncells, arma::fill::zeros); // * ncells x 1
    // ! Vetor solução   
    this->u = arma::vec(ncells, arma::fill::zeros);  // * ncells x 1
    // ! Vetor b modificado contendo correções da não ortogonalidade
    this->b_with_cd = arma::vec(ncells, arma::fill::zeros); // * ncells x 1
    
    // ! Vetores armazenando de forma pré-calculada os valores nas faces das propriedades.
    this->source = arma::vec(ncells, arma::fill::zeros);  // * nfaces
    this->gammaf = arma::vec(nfaces, arma::fill::zeros);  // * nfaces
    this->rhof = arma::vec(nfaces, arma::fill::zeros); // * nfaces
    this->uf = arma::mat(nfaces,2, arma::fill::zeros); // * nfaces
    this->pre_processing(); // * Faz o pré-processamento dos vetores acima.
}

FVMSolver::~FVMSolver(){
    //destrutor
}

/**
 * * Realiza o pré-processamento das propriedades do problema.
 */
void FVMSolver::pre_processing(){

    int nfaces = mesh->get_nedges();
    vector<Edge*> edges = mesh->get_edges();

    // ! No caso do termo fonte calcula 1 para cada célula 
    vector<Cell*> cells = mesh->get_cells();
    for(int i = 0; i < mesh->get_ncells(); i++){
        pair<double,double>& centroid = cells[i]->get_centroid();
        this->source[i] = this->sourcefunc(centroid.first, centroid.second);
    }

    for(int i = 0; i < nfaces; i++){
        // pega os nós que compartilham aquela face
        const Node* from = edges[i]->from;
        const Node* to = edges[i]->to;
        
        // pega o centro da face
        pair<double,double> midface = edges[i]->get_middle();
        double totaldist = edges[i]->get_length(); // tamanho total da aresta
        
        // distância entre meio e from
        double w1 = distance(midface.first, midface.second, to->x, to->y)/totaldist;
        
        // distância entre meio e to
        double w2 = distance(midface.first, midface.second, from->x, from->y)/totaldist;

        // * LÓGICA GERAL: aplica a função no node e interpola os valores para a face.
        this->gammaf[i] = gammafunc(to->x, to->y)*w2 + gammafunc(from->x, from->y)*w1;  
        this->rhof[i] = rhofunc(to->x, to->y)*w2 + rhofunc(from->x, from->y)*w1;  
        this->uf(i,0) = ufunc(to->x, to->y).first*w2 + ufunc(from->x, from->y).first*w1;  
        this->uf(i,1) = ufunc(to->x, to->y).second*w2 + ufunc(from->x, from->y).second*w1;  
    }

}

/**
 * * Função para descobrir o vetor n1
 * @param: p é a celula que está sendo avaliada
 * @param: n é a celula vizinha
 * @param: normal é a normal da face apontando para fora da célula
 */
pair<double,double> compute_n1(pair<double,double>& p, pair<double,double>& n, pair<double,double>& normal){
    /* vetor que une P e N */
    pair<double,double> dpn = make_pair(n.first - p.first, n.second - p.second);
    
    /* norma do vetor normal */
    double norm_normal_square = normal.first*normal.first + normal.second*normal.second; // |nf|²
    
    /* produto escalar entre normal e dpn -> nf . dpn*/
    double nfdotdpn = normal.first*dpn.first + normal.second*dpn.second; 
    
    /* n_1 = dPN * (|n_f|^2 / (n_f . dPN))*/
    double n1x = dpn.first * (norm_normal_square/nfdotdpn);
    double n1y = dpn.second * (norm_normal_square/nfdotdpn);

    return make_pair(n1x, n1y);
}

// ! Acrescenta a difusão centrada em A e b
void FVMSolver::diffusion(){
    vector<Cell*>& cells = mesh->get_cells();
    for(int i = 0; i < cells.size(); i++){ // * PARA CADA CÉLULA
        this->diffusion_of_cell(cells[i]); // ? CALCULA A DIFUSÃO DA CÉLULA ADICIONANDO na A e b
    }
}

// ! Acrescenta a convecção upwind em A e b
void FVMSolver::convection(){
    vector<Cell*>& cells = mesh->get_cells();
    for(int i = 0; i < cells.size(); i++){ // * PARA CADA CÉLULA
        this->convection_of_cell(cells[i]); // ? CALCULA A CONVECÇÃO DA CÉLULA ADICIONANDO na A e b
    }
}

void FVMSolver::convection_of_cell(Cell *c)
{
    vector<Edge *> faces = c->get_edges(); // faces/lados da célula
    vector<int> &nsigns = c->get_nsigns();
    for (int i = 0; i < faces.size(); i++)
    {
        Edge *face = faces[i];
        int nsign = nsigns[i]; // sinal corrigida da normal

        // centro da face & normal da face corrigida
        pair<double, double> &middleFace = face->get_middle();
        pair<double, double> &normal = face->get_normal();
        pair<double, double> normal_corrected = make_pair(normal.first * nsign, normal.second * nsign);

        double UdotNormal = uf(face->id, 0) * normal_corrected.first + uf(face->id, 1) * normal_corrected.second;
        double G_f = rhof[face->id] * UdotNormal * face->get_length();
        if (face->is_boundary_face())
        {
            /*=-==-==-==-==-==-==-==-= contribuição em b =-==-==-==-==-==-==-==-==-=*/
            if(G_f > 0){
                A(c->id, c->id) += G_f;
            }else{
                BoundaryLocation local = BoundaryCondition::find_location(face);
                BoundaryType bt = boundaries[local]->get_type();
                b[c->id] += - G_f * boundaries[local]->apply(middleFace.first, middleFace.second);
            }
        }
        else
        {
            /*=-==-==-==-==-==-==-==-= contribuição em A =-==-==-==-==-==-==-==-==-=*/
            if(G_f > 0){
                // phiP
                A(c->id, c->id) += G_f;
            }else{
                // phiNeighbor
                pair<int, int> &idCellsShareFace = face->get_link_face_to_cell();
                int ic1 = idCellsShareFace.first;
                int ic2 = idCellsShareFace.second;

                int nb = ic1 == c->id ? ic2 : ic1; // * pegando o vizinho da célula P
                A(c->id, nb) += G_f;
            }
        }
    }
}

/**
 * Calcula a difusão da celula acrescentando tanto em A quando em b.
 * Termo fonte também já é calculado e adicionado em b por aqui.
 */
void FVMSolver::diffusion_of_cell(Cell* c){
    vector<Edge*> faces = c->get_edges(); // faces/lados da célula
    vector<int>& nsigns = c->get_nsigns();
    pair<double,double>& centroidP = c->get_centroid(); // centroide da celula P que está sendo avaliada
    for(int i = 0 ; i < faces.size(); i++){
        Edge* face = faces[i];
        int nsign = nsigns[i]; // sinal corrigida da normal

        // centro da face e normal da face corrigida
        pair<double,double>& middleFace = face->get_middle();
        pair<double, double>& normal = face->get_normal();
        pair<double,double> normal_corrected = make_pair(normal.first * nsign, normal.second * nsign);

        if(face->is_boundary_face()){
            pair<double,double> n1 = compute_n1(centroidP, middleFace, normal_corrected);
            double normn1 = sqrt(n1.first * n1.first + n1.second * n1.second); // |n1|

            BoundaryLocation local = BoundaryCondition::find_location(face);
            BoundaryType bt = boundaries[local]->get_type();
            if(bt == DIRICHLET){
                /*=-==-==-==-==-==-==-==-= contribuição em A =-==-==-==-==-==-==-==-==-=*/
                // TODO diagonal: pensar melhor sobre essa conta.
                A(c->id, c->id) += (gammaf[face->id] * face->get_length() * normn1) / face->get_df();

                /*=-==-==-==-==-==-==-==-= contribuição em b =-==-==-==-==-==-==-==-==-=*/
                b[c->id] += (gammaf[face->id] * face->get_length() * normn1 * boundaries[local]->apply(middleFace.first, middleFace.second)) / face->get_df();
            }else if(bt == NEUMANN){
                /* Contribuição direto na b somente*/
                b(c->id) += (gammaf[face->id] * face->get_length() * boundaries[local]->apply(middleFace.first, middleFace.second));
            }
            
        }else{
            // pega celulas que compartilham aquela face.
            pair<int,int>& idCellsShareFace = face->get_link_face_to_cell();

            int ic1 = idCellsShareFace.first; 
            int ic2 = idCellsShareFace.second;
            
            int nb = ic1 == c->id ? ic2 : ic1; // * pegando o vizinho da célula P

            pair<double,double>& centroidN = this->mesh->get_cells()[nb]->get_centroid(); // centroide do vizinho
            
            pair<double,double> n1 = compute_n1(centroidP, centroidN, normal_corrected);

            double normn1 = sqrt(n1.first * n1.first + n1.second * n1.second); // |n1|

            /*=-==-==-==-==-==-==-==-= contribuição em A =-==-==-==-==-==-==-==-==-=*/
            // ! Diagonal
            A(c->id, c->id) += (gammaf[face->id] * face->get_length() * normn1) / face->get_df();
            // ! Off-diagonal
            A(c->id, nb) += - (gammaf[face->id] * face->get_length() * normn1) / face->get_df();
        }
    }

    /*=-==-==-==-==-==-==-==-= contribuição em b (TERMO FONTE) =-==-==-==-==-==-==-==-==-=*/
    b[c->id] += this->source[c->id] * c->get_area(); // S * DeltaA
}

/*
    * Uso da reconstrução com Least Squares, e interpolação da phi_f será feita 
    * usando os valores dos centroides.
    ! Least squares irá resolver Mx = y.
*/
void FVMSolver::compute_gradients(){
    vector<Cell*>& cells = mesh->get_cells();
    this->gradients = vector<pair<double, double>>(cells.size(), std::make_pair(0.0, 0.0));

    for(int i = 0; i < cells.size(); i++){ // para cada célula
        vector<Edge*>& edgesOfCell = cells[i]->get_edges();
        
        // criação da M e da y.
        arma::mat M(edgesOfCell.size(), 2);
        arma::vec y(edgesOfCell.size());
        
        for(int j = 0; j < edgesOfCell.size(); j++){ // para cada face
            int gface = edgesOfCell[j]->id;

            if(edgesOfCell[j]->is_boundary_face()){
                // Contribuição da face de contorno
                pair<double,double>& centroid = cells[i]->get_centroid();
                pair<double,double>& middleFace = edgesOfCell[j]->get_middle();
                
                double dx = middleFace.first - centroid.first;
                double dy = middleFace.second - centroid.second; 
                
                // [dx dy]
                M(j,0) = dx;
                M(j, 1) = dy;

                // [u_B - u_P]
                BoundaryLocation local = BoundaryCondition::find_location(edgesOfCell[j]);
                y[j] = boundaries[local]->apply(middleFace.first, middleFace.second) - u[i];  
            }else{
                pair<int,int>& idCellsShareFace = edgesOfCell[j]->get_link_face_to_cell();
                    
                int ic1 = idCellsShareFace.first; 
                int ic2 = idCellsShareFace.second;
                
                int nb = ic1 == i ? ic2 : ic1; //pegando o vizinho da célula P
                
                pair<double,double>& centroidP = cells[i]->get_centroid();
                pair<double,double>& centroidN = cells[nb]->get_centroid();

                double dx = centroidN.first - centroidP.first;
                double dy = centroidN.second - centroidP.second;
                
                // [dx dy]
                M(j,0) = dx;
                M(j, 1) = dy;

                // [u_N - u_P]
                y[j] = u[nb] - u[i];  
            }
        }

        // Resolve sistema sobre-determinado M x = y.
        arma::vec x = arma::solve(M, y);
        this->gradients[i] = make_pair(x[0], x[1]);
    }
}

/*
    Cômputo da difusão cruzada após obtenção dos gradientes obtidos pelo método de reconstrução.
*/
void FVMSolver::compute_cross_diffusion(){
    vector<Cell*> cells = mesh->get_cells();
    b_with_cd.zeros(); // ? zera para a próxima iteração.
    
    for(int i = 0; i < cells.size(); i++){ 
        vector<Edge*> facesOfCell = cells[i]->get_edges();
        vector<int>& nsigns = cells[i]->get_nsigns();
        
        for(int j = 0; j < facesOfCell.size(); j++){ // para cada face
            int gface = facesOfCell[j]->id;
            
            if(facesOfCell[j]->is_boundary_face()){
                //face de contorno
                pair<double,double>& centroid = cells[i]->get_centroid();
                pair<double,double>& middleFace = facesOfCell[j]->get_middle();
                pair<double, double>& normal = facesOfCell[j]->get_normal();
                
                // correção da normal
                pair<double,double> normal_corrected = make_pair(normal.first *  nsigns[j], normal.second * nsigns[j]);
                
                pair<double,double> n1 = compute_n1(centroid, middleFace, normal_corrected);

                // encontrando o n2: n1 + n2 = nf => n2 = nf - n1
                double n2x = normal_corrected.first - n1.first;
                double n2y = normal_corrected.second - n1.second;

                // grad(phi) . normal 
                // estou assumindo que o gradiente da face é o gradiente da célula.
                double graddotnormal = gradients[i].first * n2x + gradients[i].second * n2y;
                
                /*=-==-==-==-==-==-==-==-= contribuição em b =-==-==-==-==-==-==-==-==-=*/
                this->b_with_cd[i] += (gammaf[gface] * facesOfCell[j]->get_length() * graddotnormal);
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
                
                double n2x = normal_corrected.first - n1.first;
                double n2y = normal_corrected.second - n1.second;

                /*Calculo da distancia de P até center(gface) & N até center(gface)*/
                pair<double,double>& midFace = facesOfCell[j]->get_middle();
                double d1 = distance(centroidP.first, centroidP.second, midFace.first, midFace.second);
                double d2 = distance(centroidN.first, centroidN.second, midFace.first, midFace.second);

                // interpolando grad_P e grad_N para obtenção do gradiente na face.
                double gradx = (gradients[i].first*d2 + gradients[nb].first*d1)/(d1+d2);
                double grady = (gradients[i].second*d2 + gradients[nb].second*d1)/(d1+d2);

                double graddotnormal = gradx * n2x + grady * n2y;

                /*=-==-==-==-==-==-==-==-= contribuição em b =-==-==-==-==-==-==-==-==-=*/
                this->b_with_cd[i] += (gammaf[gface] * facesOfCell[j]->get_length() * graddotnormal);
            }
        }

        // !?! por fim, b_with_cd será o b (fixo) + sum(skew_P)
        this->b_with_cd[i] = this->b[i] + this->b_with_cd[i];
    }
}

void FVMSolver::solve_system(double tolerance){
    this->u = arma::spsolve(A,b);
    // double error = 1;
    // arma::vec aux;
    // int iter = 0;
    // while(error > tolerance){
    //     aux = this->b_with_cd; // obtém a cópia com valores anteriores...

    //     this->compute_gradients(); // atualiza os gradientes na reconstrução
    //     this->compute_cross_diffusion(); // usando os gradientes computa a difusão cruzada
    //     this->u = arma::spsolve(A,b_with_cd); // resolvendo usando b atualizado com a difusão cruzada. 

    //     error = arma::norm(aux - this->b_with_cd, "inf");
    //     iter += 1;
    // }

    // cout << "Convergiu em :" << iter << " iterações.\n";
    // cout << "\nSolução obtida:\n";
}


void FVMSolver::save_solution(string filename){
    
    ofstream vtk_file(filename);

    if (!vtk_file.is_open()) {
        std::cerr << "ERROR: failed to create file." << filename << "'.\n";
        exit(1);
    }

    /*Informações padrão*/
    vtk_file << "# vtk DataFile Version 3.0\n";
    vtk_file << "Malha 2D com triângulos e temperatura\n";
    vtk_file << "ASCII\n";
    vtk_file << "DATASET UNSTRUCTURED_GRID\n\n";

    vector<Node*> nodes = this->mesh->get_nodes();
    vtk_file << "POINTS " << nodes.size() << " double" << endl;
    for(int i = 0; i < nodes.size(); i++){
        vtk_file << nodes[i]->x << " " << nodes[i]->y << " " << 0 << endl;
    }
    vtk_file << endl;

    vector<Cell*> cells = mesh->get_cells();
    vtk_file << "CELLS " << cells.size() << " ";
    int sum = 0; // vai contar a quantidade de valores a serem lidos depois...
    for(int i = 0; i < cells.size(); i++)
    {
        if(cells[i]->cellType == 2)
            sum += 4;
        else if(cells[i]->cellType == 3)
            sum += 5;
    }
    vtk_file << sum << endl;
    
    for(int i = 0; i < cells.size(); i++){
        vector<Node*>& nodesOfCell = cells[i]->get_nodes();
        vtk_file << nodesOfCell.size(); 
        for(int j = 0; j < nodesOfCell.size(); j++)
            vtk_file << " " << nodesOfCell[j]->id;
        vtk_file << endl;
    }
    vtk_file << endl;

    vtk_file << "CELL_TYPES " << cells.size() << endl;
    for(int i = 0; i < cells.size(); i++){
        if(cells[i]->cellType == 2)
            vtk_file << "5" << endl;
        else if(cells[i]->cellType == 3)
            vtk_file << "9" << endl;
    }
    vtk_file << endl;
    
    vtk_file << "CELL_DATA " << cells.size() << endl;
    vtk_file << "SCALARS Temperatura double 1\n";
    vtk_file << "LOOKUP_TABLE default\n";
    for(int i = 0; i < cells.size(); i++){
        vtk_file << this->u[i] << endl;
    }

    vtk_file.close();
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