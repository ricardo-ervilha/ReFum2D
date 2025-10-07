#include "NSSolver.h"
#include "Mesh.h"
#include "Edge.h"
#include "Cell.h"
#include "BoundaryCondition.h"

NSSolver::NSSolver(Mesh *mesh, float mu, float rho){
    this->mesh = mesh;
    this->mu = mu;
    this->rho = rho;

    // ! Inicialização de todas as variáveis que NSSolver guarda.
    int ncells = this->mesh->get_ncells();
    int nfaces = this->mesh->get_nedges();

    this->A_mom = arma::sp_mat(ncells, ncells);
    this->b_mom_x = arma::vec(ncells, arma::fill::zeros);
    this->b_mom_y = arma::vec(ncells, arma::fill::zeros);

    this->A_pc = arma::sp_mat(ncells, ncells);
    this->b_pc = arma::vec(ncells, arma::fill::zeros);
    
    // velocidades corrigidas nos centroides
    this->uc = arma::vec(ncells, arma::fill::zeros);
    this->vc = arma::vec(ncells, arma::fill::zeros);
    this->pc = arma::vec(ncells, arma::fill::zeros);
    
    this->uc_aux = arma::vec(ncells, arma::fill::zeros);
    this->vc_aux = arma::vec(ncells, arma::fill::zeros);
    this->pc_aux = arma::vec(ncells, arma::fill::zeros);

    this->u_face = arma::vec(nfaces, arma::fill::zeros);
    this->v_face = arma::vec(nfaces, arma::fill::zeros);
    this->p_face = arma::vec(nfaces, arma::fill::zeros);
    
    this->mdotf = arma::vec(nfaces, arma::fill::zeros);
    this->ap = arma::vec(ncells, arma::fill::zeros);
    
    this->mdotfcorr = arma::vec(nfaces, arma::fill::zeros);
    this->pfcorr = arma::vec(nfaces, arma::fill::zeros);
    
    this->pcorr = arma::vec(ncells, arma::fill::zeros);
    this->ucorr = arma::vec(ncells, arma::fill::zeros);
    this->vcorr = arma::vec(ncells, arma::fill::zeros);
    
    // =======================================================================================

    // TODO: Preenche boundary top do lid-driven cavity flow.
    vector<Edge*> faces = mesh->get_edges();
    for(int i = 0; i < nfaces; i++){
        Edge* face = faces[i];
        if(face->from->y == 1.0 && face->to->y == 1.0){
            // Top 
            u_face[face->id] = 1.0; 
        }
    }
}   

NSSolver::~NSSolver(){
    // nada.
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
            
            if(face->is_boundary_face()){
                // * Se uma face é de contorno, então ela tem valor prescrito de velocidade => Mas, a velocidade produto escalar normal nos contornos é sempre zero, então o termo advectivo cancela => Soma difusão no aP e passsa pro outro lado difusão pro termo fonte.
                A_mom(ic, ic) = A_mom(ic,ic) + Df;
                b_mom_x[ic] = b_mom_x[ic] + u_face[idf] * Df;
                b_mom_y[ic] = b_mom_y[ic] + v_face[idf] * Df;
            }else{
                int N = get_neighbor(face->get_link_face_to_cell(), ic);
                A_mom(ic,ic) = A_mom(ic,ic) + Df + max(mf, 0.0);
                A_mom(ic,N) = -Df - max(-mf,0.0);
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
            p_face[idf] = pc[neighbor];
        }else{
            p_face[idf] = 0.5 * (pc[ic1] + pc[ic2]);
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

void NSSolver::solve_x_mom(){
    uc = arma::spsolve(A_mom, b_mom_x);
}


void NSSolver::solve_y_mom(){
    vc = arma::spsolve(A_mom, b_mom_y);
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

            double velf_x = 0.5*uc[c1] + 0.5*uc[c2];
            double velf_y = 0.5*vc[c1] + 0.5*vc[c2];
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
            
            velf_x = 0.5 * v0dp0_x/ap[c1] + 0.5 * v1dp1_x/ap[c2];
            velf_y = 0.5 * v0dp0_y/ap[c1] + 0.5 * v1dp1_y/ap[c2];
            
            double velf_p = velf_x * face->get_normal().first + velf_y * face->get_normal().second;
            
            double vdotn = vel_f_i + velf_p - (0.5*cells[c1]->get_area()/ap[c1] + 0.5*cells[c2]->get_area()/ap[c2]) * (pc[c2] - pc[c1])/face->get_df();
            
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
                
                A_pc(ic,ic) = A_pc(ic,ic) + rho * face->get_length() * (0.5 * cells[ic]->get_area()/ap[ic] + 0.5 * cells[N]->get_area()/ap[N])/face->get_df();

                A_pc(ic, N) = -rho * face->get_length() * (0.5 * cells[ic]->get_area()/ap[ic] + 0.5 * cells[N]->get_area()/ap[N])/face->get_df();
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
            pfcorr[idf] = pcorr[neighbor];
        }else{
            pfcorr[idf] = 0.5 * (pcorr[ic1] + pcorr[ic2]);
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
        if(face->is_boundary_face()) continue;
        else{
            pair<int,int> nodes_share_face = face->get_link_face_to_cell();
            int c1 = nodes_share_face.first;
            int c2 = nodes_share_face.second;

            double coeff = 0.5*cells[c1]->get_area()/ap[c1] + 0.5 * cells[c2]->get_area()/ap[c2];
            mdotfcorr[face->id] = rho*coeff*face->get_length()*(pcorr[c1]-pcorr[c2])/face->get_df();
            mdotf[face->id] = mdotf[face->id] + mdotfcorr[face->id];
        }
    }

    cout << "Convergência de u: " << arma::norm(uc-uc_aux, "inf") << endl;
    cout << "Convergência de v: " << arma::norm(vc-vc_aux, "inf") << endl;
    // !Avalia continuidade
    // double validation = 0;
    // for(int i = 0; i < cells.size(); ++i){
    //     Cell* c = cells[i];
    //     int ic = c->id;
        

    //     vector<Edge*> faces = c->get_edges();
    //     vector<int> &nsigns = c->get_nsigns();

    //     double sum_flux = 0;

    //     for(int j = 0; j < faces.size(); ++j){ // loop over faces of cell
    //         int sign = nsigns[j];
    //         Edge* face = faces[j];
    //         sum_flux = sum_flux - mdotf[face->id] * sign; 
    //     }

    //     validation = validation + fabs(sum_flux);
    // }
    // cout << "# Continuidade: " << validation << endl;
};

void NSSolver::pres_correct(double lambda_p) {
    vector<Cell*> cells = mesh->get_cells();
    
    pc_aux = pc; // * salva valores antigos
    for(int i = 0; i < cells.size(); i++){
        int ic = cells[i]->id;
        pc[ic] = pc[ic] + lambda_p * pcorr[ic];
    }
    cout << "Convergência de p: " << arma::norm(pc-pc_aux, "inf") << endl;
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