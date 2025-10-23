#include "NSSolver.h"
#include "Mesh.h"
#include "Edge.h"
#include "Cell.h"
#include "utils.h"

NSSolver::NSSolver(Mesh *mesh, float mu, float rho, vector<BoundaryCondition> bcsU, vector<BoundaryCondition> bcsV, vector<BoundaryCondition> bcsP){ 
    // + Variáveis que o usuário informa.
    this->mesh = mesh;
    this->mu = mu;
    this->rho = rho;

    // + Inicialização de todas as variáveis que NSSolver guarda.
    int ncells = this->mesh->get_ncells();
    int nfaces = this->mesh->get_nedges();
    
    // * guarda informações sobre boundary condition: {tipo e valor}.
    this->u_boundary.resize(nfaces, {NONE, 0.0});
    this->v_boundary.resize(nfaces, {NONE, 0.0});
    this->p_boundary.resize(nfaces, {NONE, 0.0});

    // matriz A de momemnto
    this->A_mom = arma::sp_mat(ncells, ncells);
    
    // matrizes b de momento
    this->b_mom_x = arma::vec(ncells, arma::fill::zeros);
    this->b_mom_y = arma::vec(ncells, arma::fill::zeros);

    // matriz A de correção de pressão
    this->A_pc = arma::sp_mat(ncells, ncells);

    // matriz b de correção de pressão
    this->b_pc = arma::vec(ncells, arma::fill::zeros);
    
    // velocidades nos centroides
    this->uc = arma::vec(ncells, arma::fill::zeros);
    this->vc = arma::vec(ncells, arma::fill::zeros);
    this->pc = arma::vec(ncells, arma::fill::zeros);

    // + Inicializando dt eu mesmo.
    this->dt = 0.0;
    // + Os vetores abaixo eu setei a condição inicial tudo zero para testar, mas o ideal é deixar isso customizável.
    this->uc_old = arma::vec(ncells, arma::fill::zeros);
    this->vc_old = arma::vec(ncells, arma::fill::zeros);
    this->pc_old = arma::vec(ncells, arma::fill::zeros);
    
    // vetores auxiliares para ajudar a validar a convergência
    this->uc_aux = arma::vec(ncells, arma::fill::zeros);
    this->vc_aux = arma::vec(ncells, arma::fill::zeros);
    this->pc_aux = arma::vec(ncells, arma::fill::zeros);

    // velocidades nas faces
    this->u_face = arma::vec(nfaces, arma::fill::zeros);
    this->v_face = arma::vec(nfaces, arma::fill::zeros);
    this->p_face = arma::vec(nfaces, arma::fill::zeros);
    
    // fluxo de massa nas faces
    this->mdotf = arma::vec(nfaces, arma::fill::zeros);
    
    // valor dos coeficientes da diagonal!
    this->ap = arma::vec(ncells, arma::fill::zeros);
    
    // correções
    this->mdotfcorr = arma::vec(nfaces, arma::fill::zeros);
    this->pfcorr = arma::vec(nfaces, arma::fill::zeros);
    this->pcorr = arma::vec(ncells, arma::fill::zeros);
    this->ucorr = arma::vec(ncells, arma::fill::zeros);
    this->vcorr = arma::vec(ncells, arma::fill::zeros);
    
    // =======================================================================================

    // chama para inicializar as condições de contorno
    bcsu = bcsU;
    bcsv = bcsV;
    bcsp = bcsP;
    this->compute_bcs();

    // chama para calcular os pesos que precisa na hora de fazer a media do valor:
    // phi_f = wf*phi_P + (1-wf)*phi_N
    this->wf = arma::vec(nfaces, arma::fill::zeros);
    this->compute_wf();
}   

NSSolver::~NSSolver(){
    // nada.
}

void NSSolver::compute_bcs(){
    
    int nfaces = mesh->get_nedges();
    vector<Edge*> faces = mesh->get_edges();

    auto apply_bc = [](auto& boundary, auto& bcs, Edge* face) {
    auto [x, y] = face->get_middle();
        for(auto& bc : bcs) {
            if(bc.get_location(x, y)) {
                boundary[face->id].first  = bc.get_type();
                boundary[face->id].second = bc.apply(x, y);
                break;
            }
        }
    };
    
    // aplica de fato a BC salvando no .first o tipo (DIRICHLET ou NEUMANN) e o valor dela, se é zero, um ou uma função que o usuário tiver passado ele aplica a função.
    for(int i = 0; i < nfaces; i++){
        Edge* face = faces[i];
        if(!face->is_boundary_face()) continue;

        apply_bc(u_boundary, bcsu, face);
        apply_bc(v_boundary, bcsv, face);
        apply_bc(p_boundary, bcsp, face);
    }

    // faz os valores das faces quando for DIRICHLET em velocidade (wall ou inlet) serem o valor prescrito e mesma coisa com pressão quando for outlet
    for(int i = 0; i < nfaces; i++){
        Edge* face = faces[i];
        int idf = face->id;
        if(face->is_boundary_face()){
            if(u_boundary[idf].first == DIRICHLET)
                u_face[idf] = u_boundary[idf].second;
                
            if(v_boundary[idf].first == DIRICHLET)
                v_face[idf] = v_boundary[idf].second;
            
            if(p_boundary[idf].first == DIRICHLET)
                p_face[idf] = p_boundary[idf].second;
        }
    }

    // Calculando fluxo de massa (Caso seja inlet irá pegar fluxo de massa diferente de zero)
    for(int i = 0; i < nfaces; i++){
        Edge* face = faces[i];
        int idf = face->id;
        pair<double,double> normal = face->get_normal();
        // ρ_f * (U . n)_f * A_f 
        // * se for condição wall => mdotf = 0
        // * se for condição inlet => mdotf != 0
        // * se for condição outlet => começa como zero, pois a velocidade na condição de contorno é igual a velocidade na celula vizinha, e a velocidade na celula vizinha é inicializa com zero, portanto mdotf = 0.
        mdotf[idf] = rho * (u_face[idf] * normal.first + v_face[idf] * normal.second) * face->get_length();
    }

    // ! IMPRIME INFORMAÇÃO NAS FACES.
    // for(int i = 0; i < nfaces; i++){
    //     cout << "u_face: " << u_face[i] << " " << " v_face: " << v_face[i] << " p_face: " << p_face[i] << " mdotf: " << mdotf[i] << endl;
    // }
    
    // for(int i = 0; i < nfaces; i++){
    //     Edge* face = faces[i];
    //     if(face->is_boundary_face())
    //         cout << "x: " << face->get_middle().first << "\ty: " << face->get_middle().second << "\t" << u_boundary[i].first << "\t" << u_boundary[i].second << endl;
    // }
    // cout << "===========================\n";

    // for(int i = 0; i < nfaces; i++){
    //     Edge* face = faces[i];
    //     if(face->is_boundary_face())
    //         cout << "x: " << face->get_middle().first << "\ty: " << face->get_middle().second << "\t" << v_boundary[i].first << "\t" << v_boundary[i].second << endl;
    // }
    // cout << "===========================\n";

    // for(int i = 0; i < nfaces; i++){
    //     Edge* face = faces[i];
    //     if(face->is_boundary_face())
    //         cout << "x: " << face->get_middle().first << "\ty: " << face->get_middle().second << "\t" << p_boundary[i].first << "\t" << p_boundary[i].second << endl;
    // }
    // cout << "===========================\n";
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
    
    // zera para recalcular a matriz.
    A_mom.zeros();
    
    /*
    * * Coeficientes de Convecção-difusão.
    */
    for(int i = 0; i < cells.size(); ++i){
        // * warm-up
        Cell* c = cells[i];
        int ic = c->id;
        vector<int> &nsigns = c->get_nsigns();
        vector<Edge*> faces = c->get_edges();
        
        // zera tudo
        b_mom_x[ic] = 0;
        b_mom_y[ic] = 0;
        
        for(int j = 0; j < faces.size(); ++j){
            // * warm-up
            Edge* face = faces[j];
            int idf = face->id;
            int nsign = nsigns[j];
            
            // Df = μ_f * A_f / δ_f 
            double Df = (mu * face->get_length()) / face->get_df();
            double mf = mdotf[idf] * nsign;            
            
            if(face->is_boundary_face()){
                
                if(u_boundary[idf].first == DIRICHLET){
                    // aP = ap + Df
                    // aN = 0
                    // S = S + Df * u_f - mf * u_f (entra a parte no termo fonte)
                    A_mom(ic, ic) = A_mom(ic,ic) + Df; // + SOMA ISSO UMA VEZ SÓ, NO PRÓXIMO NÃO É NECESSÁRIO.
                    b_mom_x[ic] = b_mom_x[ic] + u_face[idf] * Df - u_face[idf] * mf;
                } else{
                    // considera-se que u_b = u_P; v_b = v_P.
                    // com isso, o termo difusivo na discretização:  (μ_f*A_f/δ_f)*(u_b - u_P) = (μ_f*A_f/δ_f)*(u_P - u_P) = 0.
                    pair<double,double> nc = {
                        face->get_normal().first * nsign,
                        face->get_normal().second * nsign
                    };
                    // mdotf será calculado de forma alternativa, considerando a celula que faz divisa: ρ_f * A_f  * (U_c . n)
                    double mf_O = rho * face->get_length() * (nc.first * uc[ic] + nc.second * vc[ic]);
                    // difusão é zero, e termo convectivo passa pro lado do termo fonte.
                    b_mom_x[ic] = b_mom_x[ic] - mf_O * uc[ic];
                }  
                
                if(v_boundary[idf].first == DIRICHLET){
                    // aP = ap + Df
                    // aN = -Df
                    // S = S + Df * v_f - mf * v_f (convecção passou para o lado do termo fonte.)
                    b_mom_y[ic] = b_mom_y[ic] + v_face[idf] * Df - v_face[idf] * mf;
                } else{
                    // considera-se que u_b = u_P; v_b = v_P.
                    // com isso, o termo difusivo na discretização:  (μ_f*A_f/δ_f)*(v_b - v_P) = (μ_f*A_f/δ_f)*(v_P - v_P) = 0.
                    pair<double,double> nc = {
                        face->get_normal().first * nsign,
                        face->get_normal().second * nsign
                    };
                    // mdotf será calculado de forma alternativa, considerando a celula que faz divisa: ρ_f * A_f  * (U_c . n)
                    double mf_O = rho * face->get_length() * (nc.first * uc[ic] + nc.second * vc[ic]);
                    // difusão é zero, e termo convectivo passa pro lado do termo fonte.
                    b_mom_y[ic] = b_mom_y[ic] - mf_O * vc[ic];
                }
            }else{
                // interior cells
                // aP = aP + Df + max(0, mf)
                // aN = -Df - max(0, -mf)
                int N = get_neighbor(face->get_link_face_to_cell(), ic);
                
                A_mom(ic,ic) = A_mom(ic,ic) + Df + max(mf, 0.0);
                
                A_mom(ic,N) = -Df - max(-mf,0.0);
            }
        }

        // * Adição do termo transiente em aP: (ρ * V) / ∆t
        // A_mom(ic, ic) = A_mom(ic,ic) + (rho * c->get_area())/dt;

        // * Adição do termo transiente nos termos fonte: (u|v)^(n) * (ρ * V) / ∆t
        // b_mom_x[ic] = b_mom_x[ic] + uc_old[ic] * (rho * c->get_area())/dt;
        // b_mom_y[ic] = b_mom_y[ic] + vc_old[ic] * (rho * c->get_area())/dt;

        //* Após o cálculo, aplica-se a relaxação com lambda e salvar o coeficiente da diagonal para uso posterior.
        // aP <- aP / λ
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
            if(p_boundary[idf].first == NEUMANN){
                int neighbor = ic1 != -1 ? ic1 : ic2;
                // se for ∇p = 0, então assume-se que a pressão na face(pb) é igual a do vizinho p_O.
                p_face[idf] = pc[neighbor];
            }else{
                // caso contrário é dirichlet e já tem um valor prescrito pelo problema que foi preenchido quando se inicializou as condições de contorno.
            }
        }else{
            // pressão no interior é obtido via interpolação linear.
            p_face[idf] = wf[idf] * pc[ic1] + (1-wf[idf]) * pc[ic2];
        }
    }

    /**
     * * Calcula fonte do gradiente de pressão.
     */
    double regularization = (1- lambda_v);
    for(int i = 0; i < cells.size(); ++i){
        Cell* c = cells[i];
        int ic = c->id;
        vector<Edge*> faces = c->get_edges();
        vector<int> &nsigns = c->get_nsigns();
        
        for(int j = 0; j < faces.size(); ++j){
            int nsign = nsigns[j];
            Edge* face = faces[j];
            int idf = face->id;

            // b_mom_x = -Σ_f p_f * A_f * n_{x,f}
            b_mom_x[ic] = b_mom_x[ic] - p_face[idf] * face->get_normal().first * nsign * face->get_length();
            
            // b_mom_y = -Σ_f p_f * A_f * n_{y,f}
            b_mom_y[ic] = b_mom_y[ic] - p_face[idf] * face->get_normal().second * nsign * face->get_length();
        }

        // + Não esquecer, ap já está divido por λ, por isso multiplica por (1-λ) ao invés de (1-λ)/λ.
        // & Adiciona contribuição de regularização: (1-λ) * aP * u_O^n 
        b_mom_x[ic] = b_mom_x[ic] + regularization * ap[ic] * uc[ic];
        // & Adiciona contribuição de regularização: (1-λ) * aP * v_O^n 
        b_mom_y[ic] = b_mom_y[ic] + regularization * ap[ic] * vc[ic];
    }
}

void NSSolver::solve_x_mom(){
    // Solve A x = b
    uc = arma::spsolve(A_mom, b_mom_x);
}


void NSSolver::solve_y_mom(){
    // Solve A x = b 
    vc = arma::spsolve(A_mom, b_mom_y);
}

/**
 * * Realiza a interpolação de momento para obter a velocidade nas faces. (PWIM)
 */
void NSSolver::face_velocity(){
    vector<Edge*> faces = mesh->get_edges();
    vector<Cell*> cells = mesh->get_cells();
    
    for(int i = 0; i < faces.size(); ++i){
        Edge* face = faces[i];
        int idf = face->id;
        pair<int,int> nodes_share_face = face->get_link_face_to_cell();
        int c1 = nodes_share_face.first;
        int c2 = nodes_share_face.second;

        if(face->is_boundary_face()){
            if(u_boundary[idf].first == DIRICHLET) 
                continue; // Valor já está definido (wall ou inlet) e não é necessário alterá-lo.
            else{
                // As velocidades u_face e v_face quando for neumann, irei pegar do vizinho lá no momento, então aqui não precisa alterar nada.
                continue;
            }
        }
        else{
            // * Face interior.

            // interpolação linear
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
                // Σ p_f * A_f * n_x,f onde f são as faces da celula 1
                v0dp0_x = v0dp0_x + p_face[faceC1->id] * faceC1->get_length() * faceC1->get_normal().first * nsign;
                // Σ p_f * A_f * n_y,f onde f sao as faces da celula 1
                v0dp0_y = v0dp0_y + p_face[faceC1->id] * faceC1->get_length() * faceC1->get_normal().second * nsign;
            }

            double v1dp1_x = 0;
            double v1dp1_y = 0;

            vector<Edge*> facesOfc2 = cells[c2]->get_edges();
            vector<int> nsignsc2 = cells[c2]->get_nsigns();
            for(int j = 0; j < facesOfc2.size(); j++){
                // faces of cell N, contas similares ao caso de cima
                Edge* faceC2 = facesOfc2[j];
                int nsign = nsignsc2[j];
                v1dp1_x = v1dp1_x + p_face[faceC2->id] * faceC2->get_length() * faceC2->get_normal().first * nsign;
                v1dp1_y = v1dp1_y + p_face[faceC2->id] * faceC2->get_length() * faceC2->get_normal().second * nsign;
            }
            
            // continua contas do PWIM
            velf_x = wf[idf] * v0dp0_x/ap[c1] + (1-wf[idf]) * v1dp1_x/ap[c2];
            velf_y = wf[idf] * v0dp0_y/ap[c1] + (1-wf[idf]) * v1dp1_y/ap[c2];
            
            double velf_p = velf_x * face->get_normal().first + velf_y * face->get_normal().second;
            
            double vdotn = vel_f_i + velf_p - (wf[idf]*cells[c1]->get_area()/ap[c1] + (1-wf[idf])*cells[c2]->get_area()/ap[c2]) * (pc[c2] - pc[c1])/face->get_df();
            
            // no fim, tem-se o fluxo de massa da face interpolado.
            mdotf[face->id] = vdotn * rho * face->get_length();
        }
    }
}

/**
 * * Calcula a correção da pressão (p')
 */
void NSSolver::solve_pp(bool sing_matrix){
    // *Calcula termo fonte: -Σ m_f
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
            
            if(u_boundary[idf].first == NEUMANN){
                // Análogo ao caso do momento, é preciso usar as velocidades das celulas vizinhas para encontrar o mf.
                pair<double,double> normal_corrected = {
                        face->get_normal().first * nsign,
                        face->get_normal().second * nsign
                };
                double mf_outlet = rho * (uc[ic] * normal_corrected.first + vc[ic] * normal_corrected.second) * face->get_length();
                b_pc[ic] = b_pc[ic] - mf_outlet;
            }else{
                // Quando for inlet ou wall (dirichlet) => usa o valor que já tem armazenado.
                b_pc[ic] = b_pc[ic] - mdotf[idf] * nsign;
            }
        }
    }
        
   // Calcula coeficientes da A
   A_pc.zeros();
   for(int i = 0; i < cells.size(); ++i){
        Cell* c = cells[i];
        int ic = c->id;

        vector<Edge*> faces = c->get_edges();
        vector<int> &nsigns = c->get_nsigns();

        for(int j = 0; j < faces.size(); ++j){ // loop over faces of cell
            int sign = nsigns[j];
            Edge* face = faces[j];
            int idf = face->id;

            if(face->is_boundary_face()){
                if(u_boundary[idf].first == DIRICHLET){
                    // se for wall, a correção do fluxo de massa é zero, então não tem nada para fazer.
                    // se for inlet, a correção do fluxo de massa é zero também, visto que já foi especificado, então não há nada a fazer.
                    continue;
                }else{
                    // se for neumann (outlet), então...
                    // correção do fluxo de massa em geral = rho_f * Af * (V_P/aP + V_N/aN) * (p'_O - p'_N)/delta_f
                    // tomar-se então: rho_f * Af * (V_P/aP)(p'_O - p'_b)/delta_f. Mas, p'b = 0 pois a pressão é dada.
                    // => rho_f * Af * (V_P/aP) * (p'_O)/delta_f 
                    A_pc(ic, ic) = A_pc(ic,ic) + (rho * face->get_length() * (c->get_area()/ap[ic])) / face->get_df();
                }
            }
            else{

                // faces interiores
                int N = get_neighbor(face->get_link_face_to_cell(), ic);
                
                A_pc(ic,ic) = A_pc(ic,ic) + rho * face->get_length() * (wf[idf] * cells[ic]->get_area()/ap[ic] + (1-wf[idf]) * cells[N]->get_area()/ap[N])/face->get_df();

                A_pc(ic, N) = -rho * face->get_length() * (wf[idf] * cells[ic]->get_area()/ap[ic] + (1-wf[idf]) * cells[N]->get_area()/ap[N])/face->get_df();
            }
        }
    }
    
    if(sing_matrix){
        // removendo problema da matriz singular (caso do lid, no caso do backward não deveria ser.)
        A_pc.row(0).zeros();
        A_pc(0,0) = 1;
        b_pc[0] = 0;
    }

    // resolve para obter p'.
    pcorr = arma::spsolve(A_pc, b_pc);
}

/*
* * Corrige velocidades.
*/
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
            if(p_boundary[idf].first == NEUMANN){
                // se for neumann, assumo que p'_b = p'_C
                int neighbor = ic1 != -1 ? ic1 : ic2;
                pfcorr[idf] = pcorr[neighbor];
            }else{
                // caso contrário fica zero. p' = 0 (pois é dirichlet e já está certo!)
                pfcorr[idf] = 0.0;
            }
        }else{
            // interior: interpolação linear
            pfcorr[idf] = wf[idf] * pcorr[ic1] + (1-wf[idf]) * pcorr[ic2];
        }
    }
    
    // Corrigindo valores das velocidades no centro
    vector<Cell*> cells = mesh->get_cells();
    
    // + convergência (!!!) 
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
        if(face->is_boundary_face()){
            if(u_boundary[idf].first == DIRICHLET) 
                continue; // nada a fazer.
            else{
                // nada a fazer, pois não será utilizado, pois sempre pego o vizinho da célula!
            }
        }
        else{
            // interior.
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
}


void NSSolver::TransientSimple(){
    for(int i = 0; i < 1; i++){
        
        cout << "Iteração: " << i << endl;
        // // zera todos para a próxima iteração de tempo
        // uc.zeros();
        // vc.zeros();
        // pc.zeros();
        // u_face.zeros();
        // v_face.zeros();
        // p_face.zeros();
        // mdotf.zeros();
        // ap.zeros();
        // mdotfcorr.zeros();
        // pfcorr.zeros();
        // pcorr.zeros();
        // ucorr.zeros();
        // vcorr.zeros();

        // this->compute_bcs();
        
        for(int j = 0; j < 1000; j++){
            // resolve o simple 100 iterações.
            cout << "# Calculando A_mom, b_mom_x e b_mom_y\n";
            mom_links_and_sources(0.3);
            cout << "Resolvendo para encontrar uc\n";
            solve_x_mom();
            cout << "Resolvendo para encontrar uv\n";
            solve_y_mom();
            
            cout << "# Calculando velocidade nas faces {u_f e v_f}\n";
            face_velocity();
            cout << "# Calculando correção na pressão (p')\n";
            solve_pp(true); // lid chama com true, backward facing step chama com false
            cout << "# Atualiza velocidades...\n";
            uv_correct();
            cout << "# Atualiza pressão....\n";
            pres_correct(0.05);
        }
        // faz variaveis na proxima iteração começarem como as antigas.
        uc_old = uc;
        vc_old = vc;
        pc_old = pc;
        
        string fn = "../outputs/cylinder_testecase_2d_1_steady.vtk";
        export_solution(fn);
    }
}


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