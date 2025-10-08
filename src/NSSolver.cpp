#include "NSSolver.h"
#include "Mesh.h"
#include "Edge.h"
#include "Cell.h"
#include "utils.h"

NSSolver::NSSolver(Mesh *mesh, float mu, float rho, BoundaryCondition *bc1, BoundaryCondition* bc2, BoundaryCondition* bc3, BoundaryCondition *bc4, BoundaryCondition *bc5, BoundaryCondition* bc6, BoundaryCondition* bc7, BoundaryCondition *bc8, BoundaryCondition *bc9, BoundaryCondition* bc10, BoundaryCondition* bc11, BoundaryCondition *bc12){
    this->mesh = mesh;
    this->mu = mu;
    this->rho = rho;

    this->u_bcs.emplace(bc1->get_location(), bc1);
    this->u_bcs.emplace(bc2->get_location(), bc2);
    this->u_bcs.emplace(bc3->get_location(), bc3);
    this->u_bcs.emplace(bc4->get_location(), bc4);

    this->v_bcs.emplace(bc5->get_location(), bc5);
    this->v_bcs.emplace(bc6->get_location(), bc6);
    this->v_bcs.emplace(bc7->get_location(), bc7);
    this->v_bcs.emplace(bc8->get_location(), bc8);

    this->p_bcs.emplace(bc9->get_location(), bc9);
    this->p_bcs.emplace(bc10->get_location(), bc10);
    this->p_bcs.emplace(bc11->get_location(), bc11);
    this->p_bcs.emplace(bc12->get_location(), bc12);

    BoundaryCondition::apply_mins_maxs(mesh->get_xmin(), mesh->get_xmax(), mesh->get_ymin(), mesh->get_ymax());

    // ! Inicialização de todas as variáveis que NSSolver guarda.
    int ncells = this->mesh->get_ncells();
    int nfaces = this->mesh->get_nedges();

    this->A_mom = arma::sp_mat(ncells, ncells);
    this->b_mom_x = arma::vec(ncells, arma::fill::zeros);
    this->b_mom_y = arma::vec(ncells, arma::fill::zeros);

    this->A_pc = arma::sp_mat(ncells, ncells);
    this->b_pc = arma::vec(ncells, arma::fill::zeros);
    this->p_face = arma::vec(nfaces, arma::fill::zeros);
    
    // velocidades corrigidas nos centroides
    this->uc = arma::vec(ncells, arma::fill::zeros);
    this->vc = arma::vec(ncells, arma::fill::zeros);
    this->pc = arma::vec(ncells, arma::fill::zeros);
    
    this->uc_aux = arma::vec(ncells, arma::fill::zeros);
    this->vc_aux = arma::vec(ncells, arma::fill::zeros);
    this->pc_aux = arma::vec(ncells, arma::fill::zeros);
    
    this->mdotf = arma::vec(nfaces, arma::fill::zeros);
    this->ap = arma::vec(ncells, arma::fill::zeros);
    
    this->mdotfcorr = arma::vec(nfaces, arma::fill::zeros);
    this->pfcorr = arma::vec(nfaces, arma::fill::zeros);
    
    this->pcorr = arma::vec(ncells, arma::fill::zeros);
    this->ucorr = arma::vec(ncells, arma::fill::zeros);
    this->vcorr = arma::vec(ncells, arma::fill::zeros);
    
    // =======================================================================================

    this->wf = arma::vec(nfaces, arma::fill::zeros);
    this->compute_wf();
    this->compute_bc();
}   

NSSolver::~NSSolver(){
    // nada.
}

void NSSolver::compute_bc(){
    /*Trata somente U e V.*/
    vector<Edge*> faces = mesh->get_edges();
    for(int i = 0; i < faces.size(); i++){
        Edge* face = faces[i];
        pair<double,double> middleFace = face->get_middle();
        if(face->is_boundary_face()){
            // acha o local
            BoundaryLocation local = BoundaryCondition::find_location(face);
            
            // tipo e valor p/ u velocidade
            BoundaryType bt_u = u_bcs[local]->get_type();
            double bv_u = u_bcs[local]->apply(middleFace.first, middleFace.second);
            u_face.push_back({bt_u, bv_u});

            // tipo e valor p/ v velocidade
            BoundaryType bt_v = v_bcs[local]->get_type();
            double bv_v = v_bcs[local]->apply(middleFace.first, middleFace.second);
            v_face.push_back({bt_v, bv_v});

            // tipo e valor p/ pressão
            BoundaryType bt_p = p_bcs[local]->get_type();
            double bv_p = p_bcs[local]->apply(middleFace.first, middleFace.second);
            p_face_bc.push_back({bt_p, bv_p});
        }else{
            u_face.push_back({NONE, 0.0});
            v_face.push_back({NONE, 0.0});
            p_face_bc.push_back({NONE, 0.0});
        }
    }
}

/*
* * Interpolação centrada baseada na distancia dos centroides para o centro da face.
*/
void NSSolver::compute_wf(){
    vector<Cell*> cells = mesh->get_cells();
    vector<Edge*> faces = mesh->get_edges();
    for(int i = 0; i < faces.size(); i++){
        Edge* face = faces[i];
        int idf = face->id;
        if(!face->is_boundary_face()){
            pair<int,int> nodes_share_face = face->get_link_face_to_cell();
            
            int ic1 = nodes_share_face.first;
            int ic2 = nodes_share_face.second;
            
            Cell* cell1 = cells[ic1];
            Cell* cell2 = cells[ic2];
            
            double d1 = distance(cell1->get_centroid().first, cell1->get_centroid().second, face->get_middle().first, face->get_middle().second);
            double d2 = distance(cell2->get_centroid().first, cell2->get_centroid().second, face->get_middle().first, face->get_middle().second);

            wf[idf] = d2/(d1+d2);
        }
    }
}

void NSSolver::mom_links_and_sources(double lambda_v){
    vector<Cell*> cells = mesh->get_cells();
    A_mom.zeros();
    /*
    * * Convecção-difusão.
    */
    for(int i = 0; i < cells.size(); ++i){
        // * warm-up
        Cell* c = cells[i];
        int ic = c->id;
        vector<int> &nsigns = c->get_nsigns();
        vector<Edge*> faces = c->get_edges();
        
        // zera tudo
        // A_mom.row(ic).zeros();
        b_mom_x[ic] = 0;
        b_mom_y[ic] = 0;
        
        for(int j = 0; j < faces.size(); ++j){
            // * warm-up
            Edge* face = faces[j];
            int idf = face->id;
            int nsign = nsigns[j];
   
            double Df = (mu * face->get_length()) / face->get_df();
            double mf = mdotf[idf] * nsign;    
            
            pair<double, double> &middleFace = face->get_middle();
            pair<double, double> &normal = face->get_normal();
            pair<double, double> normal_corrected = make_pair(normal.first * nsign, normal.second * nsign);
            
            if(face->is_boundary_face()){
                // * Se uma face é de contorno, então ela tem valor prescrito de velocidade => Mas, a velocidade produto escalar normal nos contornos é sempre zero, então o termo advectivo cancela => Soma difusão no aP e passsa pro outro lado difusão pro termo fonte.
                A_mom(ic, ic) = A_mom(ic,ic) + Df;

                pair<double,double> n1 = compute_n1(c->get_centroid(), middleFace, normal_corrected);
                double normn1 = sqrt(n1.first * n1.first + n1.second * n1.second); // |n1|

                if(u_face[idf].first == DIRICHLET)
                    b_mom_x[ic] = b_mom_x[ic] + u_face[idf].second * Df * normn1 - u_face[idf].second * mf;
                
                if(v_face[idf].first == DIRICHLET)
                    b_mom_y[ic] = b_mom_y[ic] + v_face[idf].second * Df * normn1 - u_face[idf].second * mf;
            }else{
                int N = get_neighbor(face->get_link_face_to_cell(), ic);
                pair<double,double>& centroidN = cells[N]->get_centroid();
                pair<double,double> n1 = compute_n1(c->get_centroid(), centroidN, normal_corrected);

                double norm_n1 = sqrt(n1.first * n1.first + n1.second * n1.second); // + |n1|
                
                A_mom(ic,ic) = A_mom(ic,ic) + Df * norm_n1 + max(mf, 0.0);
                A_mom(ic,N) = -Df * norm_n1 - max(-mf,0.0);
            }
        }

        //* Após o cálculo, eu tenho que aplicar o lambda e salvar o coeficiente da diagonal para uso posterior.
        A_mom(ic, ic) = A_mom(ic,ic) / lambda_v;
        ap[ic] = A_mom(ic, ic);
    }

    /*
    * * Calcula pressão nas faces
    */
    vector<Edge*> faces = mesh->get_edges();
    for(int i = 0; i < faces.size(); i++){
        Edge* face = faces[i];
        int idf = face->id;
        pair<int,int> id_nodes_share_face = face->get_link_face_to_cell();

        int ic1 = id_nodes_share_face.first;
        int ic2 = id_nodes_share_face.second;

        if(face->is_boundary_face()){
            int neighbor = ic1 != -1 ? ic1 : ic2;
            if(p_face_bc[idf].first == NEUMANN) // SÉRIE DE TAYLOR
                p_face[idf] = pc[neighbor] + p_face_bc[idf].second * face->get_df();
        }else{
            p_face[idf] = wf[idf] * pc[ic1] + (1-wf[idf]) * pc[ic2];
        }
    }

    /**
     * * Calcula fonte do gradiente de pressão.
     */
    double regularization = ( 1- lambda_v);
    for(int i = 0; i < cells.size(); ++i){
        Cell* c = cells[i];
        int ic = c->id;
        vector<Edge*> faces = c->get_edges();
        vector<int> &nsigns = c->get_nsigns();
        
        for(int j = 0; j < faces.size(); ++j){
            int nsign = nsigns[j];
            Edge* face = faces[j];
            int idf = face->id;

            b_mom_x[ic] = b_mom_x[ic] - p_face[idf] * face->get_normal().first * nsign * face->get_length();
            
            b_mom_y[ic] = b_mom_y[ic] - p_face[idf] * face->get_normal().second * nsign * face->get_length();
        }

        // & Adiciona contribuição de regularização
        b_mom_x[ic] = b_mom_x[ic] + regularization * ap[ic] * uc[ic];
        b_mom_y[ic] = b_mom_y[ic] + regularization * ap[ic] * vc[ic];
    }
}

arma::vec NSSolver::cross_diffusion(char var, arma::mat gradients){
    vector<Cell*> cells = mesh->get_cells();

    arma::vec phi;
    vector<pair<BoundaryType, double>> phif;
    arma::vec source_cross_diffusion = arma::vec(cells.size(), arma::fill::zeros);
    if(var == 'u')
    {
        phi = uc;
        phif = u_face;
    }else if(var == 'v'){
        phi = vc;
        phif = v_face;
    }

    for(int i = 0; i < cells.size(); i++){ 
        vector<Edge*> facesOfCell = cells[i]->get_edges();
        vector<int>& nsigns = cells[i]->get_nsigns();
        pair<double,double>& centroidP = cells[i]->get_centroid();

        for(int j = 0; j < facesOfCell.size(); j++){
            Edge* e = facesOfCell[j];
            int idf = e->id;
            
            int nsign = nsigns[j];

            pair<double, double> &middleFace = e->get_middle();
            pair<double, double>& normal = e->get_normal();
            pair<double,double> normal_corrected = make_pair(normal.first * nsign, normal.second * nsign);    

            if(e->is_boundary_face()){
                 
                if(phif[idf].first == DIRICHLET){
                    pair<double,double> n1 = compute_n1(centroidP, middleFace, normal_corrected);

                    // + encontrando o n2: n1 + n2 = nf => n2 = nf - n1
                    double n2x = normal_corrected.first - n1.first;
                    double n2y = normal_corrected.second - n1.second;

                    // estou assumindo que o gradiente da face é o gradiente da celula.
                    double graddotnormal = gradients(i,0) * n2x + gradients(i,1) * n2y;

                    source_cross_diffusion[i] += (mu * e->get_length() * graddotnormal);
                }else {
                    // * não precisa tratar
                }
            }else{
                int nb = get_neighbor(e->get_link_face_to_cell(), cells[i]->id);
                pair<double,double>& centroidN = cells[nb]->get_centroid(); // centroide do vizinho

                pair<double,double> n1 = compute_n1(centroidP, centroidN, normal_corrected);

                double n2x = normal_corrected.first - n1.first;
                double n2y = normal_corrected.second - n1.second;

                /*Calculo da distancia de P até center(gface) & N até center(gface)*/
                pair<double,double>& midFace = facesOfCell[j]->get_middle();
                double d1 = distance(centroidP.first, centroidP.second, midFace.first, midFace.second);
                double d2 = distance(centroidN.first, centroidN.second, midFace.first, midFace.second);
                
                // + Interpolação do gradiente da face
                double gradx = (gradients(i,0)*d2 + gradients(nb,0)*d1)/(d1+d2);
                double grady = (gradients(i,1)*d2 + gradients(nb,1)*d1)/(d1+d2);

                double graddotnormal = gradx * n2x + grady * n2y;

                source_cross_diffusion[i] += mu * facesOfCell[j]->get_length() * graddotnormal;
            }
        }
    }

    return source_cross_diffusion;
}

void NSSolver::solve_x_mom(){
    arma::mat gradients;
    arma::vec s_cd;
    for(int i = 0; i < 1; i++){
        cout << "Gradient\n";
        gradients = gradient_reconstruction('u');
        cout << "CD\n";
        s_cd = cross_diffusion('u', gradients);
        cout << s_cd << endl;
        uc = arma::spsolve(A_mom, b_mom_x + s_cd);
    }
}


void NSSolver::solve_y_mom(){
    arma::mat gradients;
    arma::vec s_cd;
    for(int i = 0; i < 1; i++){
        cout << "Gradient\n";
        gradients = gradient_reconstruction('v');
        cout << "CD\n";
        s_cd = cross_diffusion('v', gradients);
        vc = arma::spsolve(A_mom, b_mom_y + s_cd);
    }
}


/**
 * * Realiza a interpolação de momento para obter a velocidade nas faces.
 */
void NSSolver::face_velocity(){
    vector<Edge*> faces = mesh->get_edges();
    vector<Cell*> cells = mesh->get_cells();
    
    for(int i = 0; i < faces.size(); ++i){
        Edge* face = faces[i];
        int idf = face->id;

        if(face->is_boundary_face()) continue; // IGNORE BC
        else{
            pair<int,int> nodes_share_face = face->get_link_face_to_cell();
            int c1 = nodes_share_face.first;
            int c2 = nodes_share_face.second;

            double velf_x = wf[idf]*uc[c1] + (1-wf[idf]) * uc[c2];
            double velf_y = wf[idf]*vc[c1] + (1-wf[idf]) * vc[c2];
            double vel_f_i = velf_x * face->get_normal().first + velf_y * face->get_normal().second;

            double v0dp0_x = 0;
            double v0dp0_y = 0;

            vector<Edge*> facesOfc1 = cells[c1]->get_edges();
            vector<int> nsignsc1 = cells[c1]->get_nsigns();
            for(int j = 0; j < facesOfc1.size(); j++){
                // faces of cell O
                Edge* faceC1 = facesOfc1[j];
                int nsign = nsignsc1[j];
                v0dp0_x = v0dp0_x + p_face[faceC1->id] * faceC1->get_length() * faceC1->get_normal().first * nsign;
                v0dp0_y = v0dp0_y + p_face[faceC1->id] * faceC1->get_length() * faceC1->get_normal().second * nsign;
            }

            double v1dp1_x = 0;
            double v1dp1_y = 0;

            vector<Edge*> facesOfc2 = cells[c2]->get_edges();
            vector<int> nsignsc2 = cells[c2]->get_nsigns();
            for(int j = 0; j < facesOfc2.size(); j++){
                // faces of cell N
                Edge* faceC2 = facesOfc2[j];
                int nsign = nsignsc2[j];
                v1dp1_x = v1dp1_x + p_face[faceC2->id] * faceC2->get_length() * faceC2->get_normal().first * nsign;
                v1dp1_y = v1dp1_y + p_face[faceC2->id] * faceC2->get_length() * faceC2->get_normal().second * nsign;
            }
            
            velf_x = wf[idf] * v0dp0_x/ap[c1] + (1-wf[idf]) * v1dp1_x/ap[c2];
            velf_y = wf[idf] * v0dp0_y/ap[c1] + (1-wf[idf]) * v1dp1_y/ap[c2];
            
            double velf_p = velf_x * face->get_normal().first + velf_y * face->get_normal().second;
            
            double vdotn = vel_f_i + velf_p - (wf[idf]*cells[c1]->get_area()/ap[c1] + (1-wf[idf])*cells[c2]->get_area()/ap[c2]) * (pc[c2] - pc[c1])/face->get_df();
            
            mdotf[face->id] = vdotn * rho * face->get_length();
        }
    }
}

/**
 * * Calcula a correção da pressão (p')
 */
void NSSolver::solve_pp(){

    /*Calcula termo fonte*/
    vector<Cell*> cells = mesh->get_cells();
    for(int i = 0; i < cells.size(); ++i){
        
        // * warm-up
        Cell* c = cells[i];
        int ic = c->id;
        vector<Edge*> faces = c->get_edges();
        vector<int> nsigns = c->get_nsigns();
        b_pc[ic] = 0;
        
        for(int j = 0; j < faces.size(); ++j){
            //* warm-up
            Edge* face = faces[j];
            int idf = face->id;
            int nsign = nsigns[j];

            b_pc[ic] = b_pc[ic] - mdotf[idf] * nsign;
        }
    }
        
   // Calcula coeficientes da A
   A_pc.zeros();
   for(int i = 0; i < cells.size(); ++i){
        Cell* c = cells[i];
        int ic = c->id;
        // A_pc.row(ic).zeros();

        vector<Edge*> faces = c->get_edges();
        vector<int> &nsigns = c->get_nsigns();

        for(int j = 0; j < faces.size(); ++j){ // loop over faces of cell
            int sign = nsigns[j];
            Edge* face = faces[j];
            int idf = face->id;

            if(face->is_boundary_face()) continue;
            else{
                int N = get_neighbor(face->get_link_face_to_cell(), ic);
                
                A_pc(ic,ic) = A_pc(ic,ic) + rho * face->get_length() * (wf[idf] * cells[ic]->get_area()/ap[ic] + (1-wf[idf]) * cells[N]->get_area()/ap[N])/face->get_df();

                A_pc(ic, N) = -rho * face->get_length() * (wf[idf] * cells[ic]->get_area()/ap[ic] + (1-wf[idf]) * cells[N]->get_area()/ap[N])/face->get_df();
            }
        }
    }

    A_pc.row(0).zeros();
    A_pc(0,0) = 1;
    b_pc[0] = 0;
    
    pcorr = arma::spsolve(A_pc, b_pc);
}


void NSSolver::uv_correct() {

    // Interpola pfcorr nas faces
    vector<Edge*> faces = mesh->get_edges();
    for(int i = 0; i < faces.size(); i++){
        Edge* face = faces[i];
        int idf = face->id;
        pair<int,int> id_nodes_share_face = face->get_link_face_to_cell();
        int ic1 = id_nodes_share_face.first;
        int ic2 = id_nodes_share_face.second;

        if(face->is_boundary_face()){
            int neighbor = ic1 != -1 ? ic1 : ic2;
            if(p_face_bc[idf].first == NEUMANN)
                pfcorr[idf] = pcorr[neighbor] + p_face_bc[idf].second * face->get_df();
        }else{
            pfcorr[idf] = wf[idf] * pcorr[ic1] + (1-wf[idf]) * pcorr[ic2];
        }
    }
    
    // Corrigindo valores das velocidades no centro
    vector<Cell*> cells = mesh->get_cells();
    uc_aux = uc; // * salva valor antigo
    vc_aux = vc; // * salva valor antigo
    for(int i = 0; i < cells.size(); i++){
        Cell* c = cells[i];
        int ic = c->id;
        ucorr[ic] = 0;
        vcorr[ic] = 0;
        
        vector<Edge*> faces = c->get_edges();
        vector<int> &nsigns = c->get_nsigns();
        for(int j = 0; j < faces.size(); j++){
            Edge* face = faces[j];
            int nsign = nsigns[j];

            ucorr[ic] = ucorr[ic] + pfcorr[face->id] * face->get_normal().first * nsign * face->get_length();
            vcorr[ic] = vcorr[ic] + pfcorr[face->id] * face->get_normal().second * nsign * face->get_length();
        }
        ucorr[ic] = -ucorr[ic]/ap[ic];
        vcorr[ic] = -vcorr[ic]/ap[ic];
        
        uc[ic] = uc[ic] + ucorr[ic];
        vc[ic] = vc[ic] + vcorr[ic];
    }

    // Corrigindo valores das velocidades nas faces
    for(int i = 0; i < faces.size(); i++){
        Edge* face = faces[i];
        int idf = face->id;
        if(face->is_boundary_face()) continue;
        else{
            pair<int,int> nodes_share_face = face->get_link_face_to_cell();
            int c1 = nodes_share_face.first;
            int c2 = nodes_share_face.second;

            double coeff = wf[idf]*cells[c1]->get_area()/ap[c1] + (1-wf[idf]) * cells[c2]->get_area()/ap[c2];
            mdotfcorr[idf] = rho*coeff*face->get_length()*(pcorr[c1]-pcorr[c2])/face->get_df();
            mdotf[idf] = mdotf[idf] + mdotfcorr[idf];
        }
    }

    cout << "# Convergência de u: " << arma::norm(uc-uc_aux, "inf") << endl;
    cout << "# Convergência de v: " << arma::norm(vc-vc_aux, "inf") << endl;
};

void NSSolver::pres_correct(double lambda_p) {
    vector<Cell*> cells = mesh->get_cells();
    
    pc_aux = pc; // * salva valores antigos
    for(int i = 0; i < cells.size(); i++){
        int ic = cells[i]->id;
        pc[ic] = pc[ic] + lambda_p * pcorr[ic];
    }
    cout << "# Convergência de p: " << arma::norm(pc-pc_aux, "inf") << endl;
};


void NSSolver::export_solution(string filename){

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
    vtk_file << "VECTORS velocity double\n";
    for(int i = 0; i < cells.size(); i++){
        double u = uc[i];
        double v = vc[i];
        // * SALVA U E V.
        vtk_file << u << " " << v << " 0.0\n";
    }


    vtk_file.close();
}

arma::mat NSSolver::gradient_reconstruction(char var){
    vector<Cell*>& cells = mesh->get_cells();
    arma::mat gradients_phi(cells.size(), 2);
    arma::vec phi;
    vector<pair<BoundaryType, double>> phif;
    if(var == 'u')
    {
        phi = uc;
        phif = u_face;
    }else if(var == 'v'){
        phi = vc;
        phif = v_face;
    }
    
    for(int i = 0; i < cells.size(); i++){ 

        Cell *c = cells[i];

        vector<Edge*>& edgesOfCell = c->get_edges();
        pair<double,double>& centroidP = c->get_centroid();

        // criação da Matriz M e do vetor y.
        arma::mat M(edgesOfCell.size(), 2);
        arma::vec y(edgesOfCell.size());
        
        for(int j = 0; j < edgesOfCell.size(); j++){
            Edge* e = edgesOfCell[j];
            int idf = e->id;
            
            // # meio da face
            pair<double, double> &middleFace = e->get_middle();
            
            if(e->is_boundary_face()){
            
                if(phif[idf].first == DIRICHLET){
                    double dx = middleFace.first - centroidP.first;
                    double dy = middleFace.second - centroidP.second;
                    // + [dx dy]
                    M(j,0) = dx;
                    M(j, 1) = dy;

                     // + [u_B - u_P]
                    BoundaryLocation local = BoundaryCondition::find_location(edgesOfCell[j]);
                    y[j] = phif[idf].second - phi[i];  

                }else{
                    // * No caso de neumann já tem o gradiente...
                }
            }else{
                int nb = get_neighbor(e->get_link_face_to_cell(), c->id);
                pair<double,double>& centroidN = cells[nb]->get_centroid();

                double dx = centroidN.first - centroidP.first;
                double dy = centroidN.second - centroidP.second;
                
                // + [dx dy]
                M(j,0) = dx;
                M(j, 1) = dy;

                // + [u_N - u_P]
                y[j] = phi[nb] - phi[i];  
            }
        }
        
        // + Resolve sistema sobre-determinado M x = y.
        arma::vec x = arma::solve(M, y);
        gradients_phi(i,0) = x[0];
        gradients_phi(i,1) = x[1];
    }

    return gradients_phi;
}

/**
 * * Função para descobrir o vetor n1
 * @param: p é a celula que está sendo avaliada
 * @param: n é a celula vizinha
 * @param: normal é a normal da face apontando para fora da célula
 */
pair<double,double> NSSolver::compute_n1(pair<double,double>& p, pair<double,double>& n, pair<double,double>& normal){
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