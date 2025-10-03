#include "NSSolver2.h"
#include "Mesh.h"
#include "Edge.h"
#include "Cell.h"
#include "BoundaryCondition.h"

NSSolver2::NSSolver2(Mesh *mesh, float mu, float rho, float source_x, float source_y){
    this->mesh = mesh; // Guarda um ponteiro para acessar os valores da malha.
    this->mu = mu;
    this->rho = rho;
    this->source_x = source_x;
    this->source_y = source_y;

    // ========================================================================================
    // ========================================================================================
    // ! Inicialização de todas as variáveis que NSSolver guarda.
    int ncells = this->mesh->get_ncells();
    int nfaces = this->mesh->get_nedges();

    this->u_face = arma::vec(nfaces, arma::fill::zeros);
    this->v_face = arma::vec(nfaces, arma::fill::zeros);
    
    this->uc = arma::vec(ncells, arma::fill::zeros);
    this->vc = arma::vec(ncells, arma::fill::zeros);
    this->pc = arma::vec(ncells, arma::fill::zeros);

    this->p_face = arma::vec(nfaces, arma::fill::zeros);
    this->mdotf = arma::vec(nfaces, arma::fill::zeros);
    
    this->scx = arma::vec(ncells, arma::fill::zeros);
    this->scy = arma::vec(ncells, arma::fill::zeros);
    this->ap = arma::vec(ncells, arma::fill::zeros);
    this->anb = arma::mat(ncells, 4, arma::fill::zeros);
    this->res = arma::vec(ncells, arma::fill::zeros);
    
    this->sc_p = arma::vec(ncells, arma::fill::zeros);
    this->ap_p = arma::vec(ncells, arma::fill::zeros);
    this->anb_p = arma::mat(ncells, 4, arma::fill::zeros);
    this->res_p = arma::vec(ncells, arma::fill::zeros);  
    
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
        if(fabs(face->get_middle().second - 1.0) < 1e-6){
            // Top 
            u_face[face->id] = 1.0; // u_lid
        }
    }
}   

NSSolver2::~NSSolver2(){
    // nada.
}

void NSSolver2::mom_links_and_sources() {
    vector<Cell*> cells = mesh->get_cells();
    
    for(int i = 0; i < cells.size(); i++){
        Cell* c = cells[i];
        int ic = c->id;
        ap[ic] = 0;
        scx[ic] = 0;
        scy[ic] = 0;
        
        vector<Edge*> faces = c->get_edges();
        vector<int> &nsigns = c->get_nsigns();

        for(int j = 0; j < faces.size(); j++){ // loop over faces of cell
            int nsign = nsigns[j];
            Edge* face = faces[j];
            
            double mf = mdotf[face->id] * nsign;

            if(face->is_boundary_face()){
                ap[ic] = ap[ic] + mu * face->get_length() / face->get_df();
                anb(ic, j) = 0;
                scx[ic] = scx[ic] + u_face[face->id] * mu * face->get_length()/face->get_df() - mf * u_face[face->id];
                scy[ic] = scy[ic] + v_face[face->id] * mu * face->get_length()/face->get_df() - mf * v_face[face->id];
            }else{
                // interior
                ap[ic] = ap[ic] + mu * face->get_length() / face->get_df() + max(0.0, mf); 
                anb(ic, j) = -mu * face->get_length() / face->get_df() + max(0.0, -mf);
            }
        }
    }

    // calculate p_face
    vector<Edge*> faces = mesh->get_edges();
    for(int i = 0; i < faces.size(); i++){
        Edge* face = faces[i];
        pair<int,int> id_nodes_share_face = face->get_link_face_to_cell();

        if(face->is_boundary_face()){
            int neighbor = id_nodes_share_face.first != -1 ? id_nodes_share_face.first : id_nodes_share_face.second;
            p_face[face->id] = pc[neighbor];
        }else{
            p_face[face->id] = 0.5 * (pc[id_nodes_share_face.first] + pc[id_nodes_share_face.second]);
        }
    }

    // cout << "* maior pf: " << arma::max(p_face) << endl;
    // cout << "* menor pface: " << arma::min(p_face) << endl;

    // calculation of sources
    for(int i = 0; i < cells.size(); ++i){
        Cell* c = cells[i];
        int ic = c->id;
        vector<Edge*> faces = c->get_edges();
        vector<int> &nsigns = c->get_nsigns();
        for(int j = 0; j < faces.size(); ++j){ // loop over faces of cell
            int nsign = nsigns[j];
            Edge* face = faces[j];

            scx[ic] = scx[ic] - p_face[face->id] * face->get_normal().first * nsign * face->get_length();
            scy[ic] = scy[ic] - p_face[face->id] * face->get_normal().second * nsign * face->get_length();
        }
    }

    // cout << ap << endl;
    // cout << anb << endl;
};

void NSSolver2::solve_x_mom(int iter_mom, double rin_uv, double tol_inner) {
    vector<Cell*> cells = mesh->get_cells();
    double res2o;
    double res_u;
    double res2;
    for(int iter = 0; iter < iter_mom; iter++){
        double sumr = 0;

        for(int i = 0; i < cells.size(); ++i){
            Cell* c = cells[i];
            int ic = c->id;
            vector<Edge*> faces = c->get_edges();
            vector<int> &nsigns = c->get_nsigns();

            double sumf = 0;
            for(int j = 0; j < faces.size(); ++j){ // loop over faces of cell
                Edge* face = faces[j];

                if(!face->is_boundary_face()){
                    int N = get_neighbor(face->get_link_face_to_cell(), ic);
                    sumf = sumf + anb(ic, j) * uc[N];
                }
            }
            res[ic] = scx[ic] - ap[ic] * uc[ic] - sumf;
            sumr = sumr + res[ic]*res[ic];
        }

        res2 = sqrt(max(0.0, sumr));
        res2o;
        res_u;
        if(iter == 0){
            res2o = res2;
            res_u = res2;
        }

        for(int i = 0; i < cells.size(); ++i){
            Cell* c = cells[i];
            int ic = c->id;
            double utild = (res[ic])/(ap[ic] * (1 + rin_uv));
            uc[ic] = uc[ic] + utild;
        }

        if(res2/res2o < tol_inner) break;
    }

    cout << "Maior: " << arma::max(uc) << endl;
    cout << "Menor: " << arma::min(uc) << endl;
}

void NSSolver2::solve_y_mom(int iter_mom, double rin_uv, double tol_inner) {
    vector<Cell*> cells = mesh->get_cells();
    double res2o;
    double res_u;
    double res2;
    for(int iter = 0; iter < iter_mom; iter++){
        double sumr = 0;

        for(int i = 0; i < cells.size(); ++i){
            Cell* c = cells[i];
            int ic = c->id;
            vector<Edge*> faces = c->get_edges();
            vector<int> &nsigns = c->get_nsigns();

            double sumf = 0;
            for(int j = 0; j < faces.size(); ++j){ // loop over faces of cell
                Edge* face = faces[j];

                if(!face->is_boundary_face()){
                    int N = get_neighbor(face->get_link_face_to_cell(), ic);
                    sumf = sumf + anb(ic, j) * vc[N];
                }
            }
            res[ic] = scy[ic] - ap[ic] * vc[ic] - sumf;
            sumr = sumr + res[ic]*res[ic];
        }

        res2 = sqrt(max(0.0, sumr));
        res2o;
        res_u;
        if(iter == 0){
            res2o = res2;
            res_u = res2;
        }

        for(int i = 0; i < cells.size(); ++i){
            Cell* c = cells[i];
            int ic = c->id;
            double vtild = (res[ic])/(ap[ic] * (1 + rin_uv));
            vc[ic] = vc[ic] + vtild;
        }

        if(res2/res2o < tol_inner) break;
    }

    cout << "Maior: " << arma::max(vc) << endl;
    cout << "Menor: " << arma::min(vc) << endl;
};

void NSSolver2::face_velocity() {
    vector<Edge*> faces = mesh->get_edges();
    vector<Cell*> cells = mesh->get_cells();
    for(int i = 0; i < faces.size(); ++i){
        Edge* face = faces[i];

        if(face->is_boundary_face()) continue;
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

void NSSolver2::solve_pp(int iter_pp, double tol_inner) {
    pcorr.zeros();
    
    // calculate source
    double res2o;
    double res_pp;
    double res2;
    vector<Cell*> cells = mesh->get_cells();
    for(int i = 0; i < cells.size(); ++i){
        Cell* c = cells[i];
        int ic = c->id;
        sc_p[ic] = 0;

        vector<Edge*> faces = c->get_edges();
        vector<int> &nsigns = c->get_nsigns();

        for(int j = 0; j < faces.size(); ++j){ // loop over faces of cell
            int sign = nsigns[j];
            Edge* face = faces[j];
            sc_p[ic] = sc_p[ic] - mdotf[face->id] * sign; 
        }
    }

    // calculate links
    for(int i = 0; i < cells.size(); ++i){
        Cell* c = cells[i];
        int ic = c->id;
        ap_p[ic] = 0;

        vector<Edge*> faces = c->get_edges();
        vector<int> &nsigns = c->get_nsigns();

        for(int j = 0; j < faces.size(); ++j){ // loop over faces of cell
            int sign = nsigns[j];
            Edge* face = faces[j];

            anb_p(ic, j) = 0;

            if(face->is_boundary_face()) continue;
            else{
                int N = get_neighbor(face->get_link_face_to_cell(), ic);
                ap_p[ic] = ap_p[ic] + rho * face->get_length() * (0.5 * cells[ic]->get_area()/ap[ic] + 0.5 * cells[N]->get_area()/ap[N])/face->get_df();

                anb_p(ic, j) = - rho * face->get_length() * (0.5 * cells[ic]->get_area()/ap[ic] + 0.5 * cells[N]->get_area()/ap[N])/face->get_df();
            }
        }
    }

    // solve
    for(int iter = 0; iter < iter_pp; iter++){

        for(int i = 0; i < cells.size(); ++i){
            Cell* c = cells[i];
            int ic = c->id;
            vector<Edge*> faces = c->get_edges();
            vector<int> &nsigns = c->get_nsigns();
            double sumf = 0;
            for(int j = 0; j < faces.size(); j++){
                Edge* face = faces[j];

                if(!face->is_boundary_face()){
                    int N = get_neighbor(face->get_link_face_to_cell(), ic);
                    sumf = sumf + anb_p(ic, j) * pcorr[N];
                }
            } 
            pcorr[ic] = (sc_p[ic] - sumf)/ap_p[ic];
        }

        // residuals
        res_p.zeros();
        double sumr = 0;
        for(int i = 0; i < cells.size(); ++i){
            Cell* c = cells[i];
            int ic = c->id;
            vector<Edge*> faces = c->get_edges();
            vector<int> &nsigns = c->get_nsigns();
            double sumf = 0;
            for(int j = 0; j < faces.size(); j++){
                Edge* face = faces[j];

                if(!face->is_boundary_face()){
                    int N = get_neighbor(face->get_link_face_to_cell(), ic);
                    sumf = sumf + anb_p(ic, j) * pcorr[N];
                }
            } 
            res_p[ic] = sc_p[ic] - ap_p[ic] * pcorr[ic] - sumf;
            sumr = sumr + res_p[ic] * res_p[ic];
        }

        res2 = sqrt(max(0.0, sumr));
        
        if(iter == 0){
            res2o = res2;
            res_pp = res2;
        }

        if(res2/res2o < tol_inner) break;
    }
};

void NSSolver2::uv_correct(double relax_uv) {

    // calculcate face values of pressure correction
    vector<Edge*> faces = mesh->get_edges();
    for(int i = 0; i < faces.size(); i++){
        Edge* face = faces[i];
        pair<int,int> id_nodes_share_face = face->get_link_face_to_cell();

        if(face->is_boundary_face()){
            int neighbor = id_nodes_share_face.first != -1 ? id_nodes_share_face.first : id_nodes_share_face.second;
            pfcorr[face->id] = pcorr[neighbor];
        }else{
            pfcorr[face->id] = 0.5 * (pcorr[id_nodes_share_face.first] + pcorr[id_nodes_share_face.second]);
        }
    }
    
    // cout << "pcorr: " << arma::norm(pcorr, "inf") << endl;
    // cout << "Pfcorr norm: " << arma::norm(pfcorr, "inf")<< endl;

    // correct values of centers
    vector<Cell*> cells = mesh->get_cells();
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
        uc[ic] = uc[ic] + relax_uv*ucorr[ic];
        vc[ic] = vc[ic] + relax_uv*vcorr[ic];
    }

    cout << "#Maior u: " << arma::max(uc) << endl;
    cout << "#Menor u: " << arma::min(uc) << endl;

    cout << "#Maior v: " << arma::max(vc) << endl;
    cout << "#Menor v: " << arma::min(vc) << endl;
    
    // cout << "ucorr norm: " << arma::norm(ucorr, "inf")<< endl;
    // cout << "vcorr norm: " << arma::norm(vcorr, "inf")<< endl;

    for(int i = 0; i < faces.size(); i++){
        Edge* face = faces[i];
        if(face->is_boundary_face()) continue;
        else{
            pair<int,int> nodes_share_face = face->get_link_face_to_cell();
            int c1 = nodes_share_face.first;
            int c2 = nodes_share_face.second;

            double coeff = 0.5*cells[c1]->get_area()/ap[c1] + 0.5 * cells[c2]->get_area()/ap[c2];
            mdotfcorr[face->id] = rho*coeff*face->get_length()*(pcorr[c1]-pcorr[c2])/face->get_df();
            mdotf[face->id] = mdotf[face->id] + relax_uv*mdotfcorr[face->id];
        }
    }

    // cout << "mdotcorr: " << arma::norm(mdotfcorr, "inf") << endl;

    double validation = 0;
    for(int i = 0; i < cells.size(); ++i){
        Cell* c = cells[i];
        int ic = c->id;
        

        vector<Edge*> faces = c->get_edges();
        vector<int> &nsigns = c->get_nsigns();

        double sum_flux = 0;

        for(int j = 0; j < faces.size(); ++j){ // loop over faces of cell
            int sign = nsigns[j];
            Edge* face = faces[j];
            sum_flux = sum_flux - mdotf[face->id] * sign; 
        }

        validation = validation + fabs(sum_flux);
    }
    cout << "Continuidade: " << validation << endl;
};

void NSSolver2::pres_correct(double relax_p) {
    vector<Cell*> cells = mesh->get_cells();
    for(int i = 0; i < cells.size(); i++){
        int ic = cells[i]->id;
        pc[ic] = pc[ic] + relax_p * pcorr[ic];
    }
    // cout << "#Maior pc: " << arma::max(pc) << endl;
    // cout << "#Maior pc: " << arma::min(pc) << endl;
    cout << "==============================================================================\n";
};

void NSSolver2::export_solution(string filename){

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
        vtk_file << u << " " << v << " 0.0\n";
    }


    vtk_file.close();
}