#include "ReFumSolver.h"
#include "Mesh.h"
#include "Edge.h"
#include "Cell.h"
#include "utils.h"

ReFumSolver::ReFumSolver(Mesh *mesh, float mu, float rho, vector<BoundaryCondition> bcsU, vector<BoundaryCondition> bcsV, vector<BoundaryCondition> bcsP, SolverType solver){ 
    
    // + Variáveis que o usuário informa.
    this->mesh = mesh;
    this->mu = mu;
    this->rho = rho;

    this->solver = solver;

    // + Inicialização de todas as variáveis que ReFumSolver guarda.
    int ncells = this->mesh->get_ncells();
    int nfaces = this->mesh->get_nedges();

    A.resize(ncells, ncells);
    triplets.reserve(6*ncells);
    diags.resize(ncells);
    
    // * guarda informações sobre boundary condition em um formato mais util para o ReFum: {tipo e valor}.
    this->u_boundary.resize(nfaces, {NONE, 0.0});
    this->v_boundary.resize(nfaces, {NONE, 0.0});
    this->p_boundary.resize(nfaces, {NONE, 0.0});

    this->gradients.resize(ncells, 2); // um gradiente reconstruido para cada celula
    
    // matrizes b de momento
    this->b_mom_x.resize(ncells);
    this->b_mom_x.setZero();
    this->b_mom_y.resize(ncells);
    this->b_mom_y.setZero();

    // matriz b de correção de pressão
    this->b_pc.resize(ncells);
    this->b_pc.setZero();
    
    // velocidades nos centroides
    this->uc.resize(ncells);
    this->uc.setZero();
    this->vc.resize(ncells);
    this->vc.setZero();
    this->pc.resize(ncells);
    this->pc.setZero();

    // + Define esses vetores como zeros. 
    this->uc_old.resize(ncells);
    this->uc_old.setZero();
    this->vc_old.resize(ncells);
    this->vc_old.setZero();
    this->pc_old.resize(ncells);
    this->pc_old.setZero();

    // vetores auxiliares para ajudar a validar a convergência
    this->uc_aux.resize(ncells);
    this->uc_aux.setZero();
    this->vc_aux.resize(ncells);
    this->vc_aux.setZero();
    this->pc_aux.resize(ncells);
    this->pc_aux.setZero();

    // velocidades nas faces
    this->u_face.resize(nfaces);
    this->u_face.setZero();
    this->v_face.resize(nfaces);
    this->v_face.setZero();
    this->p_face.resize(nfaces);
    this->p_face.setZero();
    
    // fluxo de massa nas faces
    this->mdotf.resize(nfaces);
    this->mdotf.setZero();
    
    // valor dos coeficientes da diagonal!
    this->ap.resize(ncells);
    this->ap.setZero();
    
    // correções
    this->mdotfcorr.resize(nfaces);
    this->mdotfcorr.setZero();
    this->pfcorr.resize(nfaces);
    this->pfcorr.setZero();
    this->pcorr.resize(ncells);
    this->pcorr.setZero();
    this->ucorr.resize(ncells);
    this->ucorr.setZero();
    this->vcorr.resize(ncells);
    this->vcorr.setZero();

    // * PEGANDO ARRAYS DO MESH (pré-processamento)
    // ===== CELLS =====
    cnsigns    = mesh->getCellNsign();
    careas     = mesh->getCellArea();
    ccentroids = mesh->getCellCentroid();
    idFacesFromCell = mesh->getFacesFromCell();

    // ===== FACES =====
    flengths   = mesh->getFaceLength();
    fdfs       = mesh->getFaceDf();
    fmiddles  = mesh->getFaceMiddle();
    fnormals  = mesh->getFaceNormal();
    flftcs  = mesh->getFaceLftc();
    fboundaryfaces = mesh->getIsBoundaryFace();

    // =======================================================================================
    // chama para inicializar as condições de contorno
    bcsu = bcsU;
    bcsv = bcsV;
    bcsp = bcsP;
    this->compute_bcs_first();
    
    // chama para calcular os pesos que precisa na hora de fazer a media do valor:
    // phi_f = wf*phi_P + (1-wf)*phi_N
    this->wf.resize(nfaces);
    this->wf.setZero();
    this->compute_wf();

    this->init_gradients();
}   

ReFumSolver::~ReFumSolver(){
    // nada.
}

// *=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=

/**
 * * Define u_0, v_0 e p_0.
 */
void ReFumSolver::set_initial_condition(function<double(double, double)> u_func, function<double(double, double)> v_func, function<double(double, double)> p_func){
    int ncells = mesh->get_ncells();
    for(int ic = 0; ic < ncells; ++ic){
        double xc = ccentroids[ic].first, yc = ccentroids[ic].second;
        
        uc_old[ic] = u_func(xc, yc);
        vc_old[ic] = v_func(xc, yc);
        pc_old[ic] = p_func(xc, yc);
    }
}

/*
    * computa os pesos de interpolação.
*/
void ReFumSolver::compute_wf(){
    vector<Cell*> cells = mesh->get_cells();
    vector<Edge*> faces = mesh->get_edges();
    for(int idf = 0; idf < faces.size(); idf++){
        
        if(!fboundaryfaces[idf]){
            pair<int,int> nodes_share_face = flftcs[idf];
            int ic1 = nodes_share_face.first;
            int ic2 = nodes_share_face.second;
            double ic1x = ccentroids[ic1].first, ic1y = ccentroids[ic1].second;
            double ic2x = ccentroids[ic2].first, ic2y = ccentroids[ic2].second;

            double fx = fmiddles[idf].first, fy = fmiddles[idf].second;
            
            // distance from ic1 to face middle
            double d1 = distance(ic1x, ic1y, fx, fy);
            // distance from ic2 to face middle
            double d2 = distance(ic2x, ic2y, fx, fy);

            wf[idf] = d2/(d1+d2);
        }
    }
}

void ReFumSolver::compute_bcs_first(){
    // + código único
    int nfaces = mesh->get_nedges();
    vector<Edge*> faces = mesh->get_edges();

    auto apply_bc = [](auto& boundary, auto& bcs, Edge* face) {
        for(auto& bc : bcs) {
            
            if(bc.get_location() == face->get_physicalgroup()) {
                boundary[face->id].first  = bc.get_type();
                boundary[face->id].second = bc.apply(face->get_middle().first, face->get_middle().second);
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
}

void ReFumSolver::compute_bcs_repeat(){
    // + código que pode ser chamado várias vezes.
    int nfaces = mesh->get_nedges();
    vector<Edge*> faces = mesh->get_edges();

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
    for(int idf = 0; idf < nfaces; idf++){
        pair<double,double> normal = fnormals[idf];
        // ρ_f * (U . n)_f * A_f 
        // * se for condição wall => mdotf = 0
        // * se for condição inlet => mdotf != 0
        // * se for condição outlet => começa como zero, pois a velocidade na condição de contorno é igual a velocidade na celula vizinha, e a velocidade na celula vizinha é inicializa com zero, portanto mdotf = 0.
        mdotf[idf] = rho * (u_face[idf] * normal.first + v_face[idf] * normal.second) * flengths[idf];
    }
}

// *=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=

void ReFumSolver::init_gradients()
{
    int ncells = mesh->get_ncells();
    grad_weights.resize(ncells);
    for (int ic = 0; ic < ncells; ic++)
    {
        vector<int> fids = idFacesFromCell[ic];
        pair<double, double> centroidP = ccentroids[ic];

        Eigen::MatrixXd M(fids.size(), 2);

        for (int j = 0; j < fids.size(); j++)
        {
            int idf = fids[j];

            // # meio da face
            pair<double, double> middleFace = fmiddles[idf];

            if (fboundaryfaces[idf])
            {

                if (u_boundary[idf].first == DIRICHLET)
                {

                    double dx = middleFace.first - centroidP.first;
                    double dy = middleFace.second - centroidP.second;
                    // + [dx dy]
                    M(j, 0) = dx;
                    M(j, 1) = dy;
                }
                else if (u_boundary[idf].first == DIRICHLET)
                {
                    continue;
                    // * No caso de neumann já tem o gradiente...
                }
            }
            else
            {
                int nb = get_neighbor(flftcs[idf], ic);
                pair<double, double> centroidN = ccentroids[nb];

                double dx = centroidN.first - centroidP.first;
                double dy = centroidN.second - centroidP.second;

                // + [dx dy]
                M(j, 0) = dx;
                M(j, 1) = dy;
            }
        }
        grad_weights[ic] = (M.transpose() * M).inverse() * M.transpose();
    }
}

void ReFumSolver::reconstruct_gradients(Eigen::VectorXd var_center, Eigen::VectorXd var_face, vector<pair<BoundaryType, double>> boundary){
    int ncells = mesh->get_ncells();
    
    for(int ic = 0; ic < ncells; ic++){ 

        vector<int> fids = idFacesFromCell[ic];

        Eigen::VectorXd y(fids.size());
        
        for(int j = 0; j < fids.size(); j++){
            int idf = fids[j];
            
            if(fboundaryfaces[idf]){
                
                if(boundary[idf].first == DIRICHLET){

                     // + [u_B - u_P]
                    y[j] = var_face[idf] - var_center[ic];  

                }else if(boundary[idf].first == DIRICHLET){
                    continue;
                    // * No caso de neumann já tem o gradiente...
                }
            }else{
                int nb = get_neighbor(flftcs[idf], ic);

                // + [u_N - u_P]
                y[j] = var_center[nb] - var_center[ic];  
            }
        }
        
        // + Resolve sistema sobre-determinado M x = y.
        Eigen::VectorXd x = grad_weights[ic] * y;
        this->gradients(ic,0) = x[0];
        this->gradients(ic,1) = x[1];
    }
}

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

Eigen::VectorXd ReFumSolver::cross_diffusion(Eigen::VectorXd b, vector<pair<BoundaryType, double>> boundary){


    int ncells = mesh->get_ncells();
    for(int ic = 0; ic < ncells; ic++){ 
        
        vector<int> fids = idFacesFromCell[ic];
        vector<int> nsigns = cnsigns[ic];
        pair<double,double> centroidP = ccentroids[ic];

        for(int j = 0; j < fids.size(); j++){
            int idf = fids[j];
            int nsign = nsigns[j];

            pair<double, double> middleFace = fmiddles[idf];
            pair<double, double> normal = fnormals[idf];
            pair<double,double> normal_corrected = make_pair(normal.first * nsign, normal.second * nsign);    

            if(fboundaryfaces[idf]){
                
                if(boundary[idf].first == DIRICHLET){
                    pair<double,double> n1 = compute_n1(centroidP, middleFace, normal_corrected);

                    // + encontrando o n2: n1 + n2 = nf => n2 = nf - n1
                    double n2x = normal_corrected.first - n1.first;
                    double n2y = normal_corrected.second - n1.second;

                    // estou assumindo que o gradiente da face é o gradiente da celula.
                    double graddotnormal = gradients(ic,0) * n2x + gradients(ic,1) * n2y;

                    b[ic] += (mu * flengths[idf] * graddotnormal);
                }else if(boundary[idf].first == DIRICHLET){
                    continue;
                    // * não precisa tratar
                }
            }else{
                int nb = get_neighbor(flftcs[idf], ic);
                pair<double,double>& centroidN = ccentroids[nb]; // centroide do vizinho

                pair<double,double> n1 = compute_n1(centroidP, centroidN, normal_corrected);

                double n2x = normal_corrected.first - n1.first;
                double n2y = normal_corrected.second - n1.second;

                /*Calculo da distancia de P até center(gface) & N até center(gface)*/
                pair<double,double>& midFace = fmiddles[idf];
                double d1 = distance(centroidP.first, centroidP.second, midFace.first, midFace.second);
                double d2 = distance(centroidN.first, centroidN.second, midFace.first, midFace.second);
                
                // + Interpolação do gradiente da face
                double gradx = (gradients(ic,0)*d2 + gradients(nb,0)*d1)/(d1+d2);
                double grady = (gradients(ic,1)*d2 + gradients(nb,1)*d1)/(d1+d2);

                double graddotnormal = gradx * n2x + grady * n2y;

                b[ic] += (mu * flengths[idf] * graddotnormal);
            }
        }
    }

    return b;
}

void ReFumSolver::mom_links_and_sources(double lambda_v){
    int ncells = mesh->get_ncells();
    int nfaces = mesh->get_nedges();

    // zera para recalcular a matriz.
    triplets.clear();
    diags.setZero();
    
    /*
    * * Coeficientes de Convecção-difusão.
    */

    for(int ic = 0; ic < ncells; ++ic){
        // * warm-up
        vector<int> nsigns = cnsigns[ic];
        vector<int> fids = idFacesFromCell[ic];
        
        // zera valores anteriores
        b_mom_x[ic] = 0;
        b_mom_y[ic] = 0;
        
        for(int j = 0; j < fids.size(); ++j){
            // * warm-up
            int idf = fids[j];
            int nsign = nsigns[j];
            
            // Df = μ_f * A_f / δ_f 
            double Df = (mu * flengths[idf]) / fdfs[idf];
            double mf = mdotf[idf] * nsign;            
            
            
            pair<double,double> centroidP = ccentroids[ic];
            pair<double, double> middleFace = fmiddles[idf];
            pair<double, double> &normal = fnormals[idf];
            pair<double, double> normal_corrected = make_pair(normal.first * nsign, normal.second * nsign);
            
            if(fboundaryfaces[idf]){
                if(u_boundary[idf].first == DIRICHLET){
                    // aP = ap + Df
                    // aN = 0
                    // S = S + Df * u_f - mf * u_f (entra a parte no termo fonte)
                    // A_mom(ic, ic) = A_mom(ic,ic) + Df; // + SOMA ISSO UMA VEZ SÓ, NO PRÓXIMO NÃO É NECESSÁRIO.
                    
                    // * Quando é dirichlet N1 ligará centro de P e centro da face.
                    
                    pair<double,double> n1 = compute_n1(centroidP, middleFace, normal_corrected);
                    double normn1 = sqrt(n1.first * n1.first + n1.second * n1.second); // |n1|

                    diags[ic] = diags[ic] + Df * normn1;
                    b_mom_x[ic] = b_mom_x[ic] + u_face[idf] * Df * normn1 - u_face[idf] * mf;
                } else{
                    // considera-se que u_b = u_P; v_b = v_P.
                    // com isso, o termo difusivo na discretização:  (μ_f*A_f/δ_f)*(u_b - u_P) = (μ_f*A_f/δ_f)*(u_P - u_P) = 0.
                    pair<double,double> nc = fnormals[idf];
                    
                    // mdotf será calculado de forma alternativa, considerando a celula que faz divisa: ρ_f * A_f  * (U_c . n)
                    double mf_O = rho * flengths[idf] * (nc.first * uc[ic] + nc.second * vc[ic]);
                    // difusão é zero, e termo convectivo passa pro lado do termo fonte.
                    b_mom_x[ic] = b_mom_x[ic] - mf_O * uc[ic];
                }  
                
                if(v_boundary[idf].first == DIRICHLET){
                    // aP = ap + Df
                    // aN = -Df
                    // S = S + Df * v_f - mf * v_f (convecção passou para o lado do termo fonte.)
                    pair<double,double> n1 = compute_n1(centroidP, middleFace, normal_corrected);
                    double normn1 = sqrt(n1.first * n1.first + n1.second * n1.second); // |n1|

                    b_mom_y[ic] = b_mom_y[ic] + v_face[idf] * Df * normn1 - v_face[idf] * mf;
                } else{
                    // considera-se que u_b = u_P; v_b = v_P.
                    // com isso, o termo difusivo na discretização:  (μ_f*A_f/δ_f)*(v_b - v_P) = (μ_f*A_f/δ_f)*(v_P - v_P) = 0.
                    pair<double,double> nc = fnormals[idf];
                    
                    // mdotf será calculado de forma alternativa, considerando a celula que faz divisa: ρ_f * A_f  * (U_c . n)
                    double mf_O = rho * flengths[idf] * (nc.first * uc[ic] + nc.second * vc[ic]);
                    // difusão é zero, e termo convectivo passa pro lado do termo fonte.
                    b_mom_y[ic] = b_mom_y[ic] - mf_O * vc[ic];
                }
            }else{
                // interior cells
                // aP = aP + Df + max(0, mf)
                // aN = -Df - max(0, -mf)
                int N = get_neighbor(flftcs[idf], ic);

                pair<double,double> centroidN = ccentroids[N]; // centroide do vizinho
            
                // + Vetor n1 unitário que liga o centroide de P com N.
                pair<double,double> n1 = compute_n1(centroidP, centroidN, normal_corrected);
                double normn1 = sqrt(n1.first * n1.first + n1.second * n1.second); // + |n1|
                
                // A_mom(ic,ic) = A_mom(ic,ic) + Df + max(mf, 0.0);
                diags[ic] = diags[ic] + Df * normn1 + max(mf, 0.0);
                
                // A_mom(ic,N) = -Df - max(-mf,0.0);
                triplets.push_back(Eigen::Triplet<double>(ic, N, -Df*normn1 - max(-mf,0.0)));
            }
        }

        if(solver == Transient){
            // * Adição do termo transiente em aP: (ρ * V) / ∆t
            // A_mom(ic, ic) = A_mom(ic,ic) + (rho * careas[ic])/dt;
            diags[ic] = diags[ic] + (rho * careas[ic])/dt;

            // * Adição do termo transiente nos termos fonte: (u|v)^(n) * (ρ * V) / ∆t
            b_mom_x[ic] = b_mom_x[ic] + uc_old[ic] * (rho * careas[ic])/dt;
            b_mom_y[ic] = b_mom_y[ic] + vc_old[ic] * (rho * careas[ic])/dt;
        }
            

        //* Após o cálculo, aplica-se a relaxação com lambda e salvar o coeficiente da diagonal para uso posterior.
        // aP <- aP / λ
        // A_mom(ic, ic) = A_mom(ic,ic) / lambda_v;
        diags[ic] = diags[ic] / lambda_v;
        
        ap[ic] = diags[ic];

        triplets.push_back(Eigen::Triplet<double>(ic, ic, diags[ic]));
    }

    /*
    * * Calcula pressão nas faces
    */
    for(int idf = 0; idf < nfaces; ++idf){
    
        pair<int,int> id_nodes_share_face = flftcs[idf];

        int ic1 = id_nodes_share_face.first;
        int ic2 = id_nodes_share_face.second;

        if(fboundaryfaces[idf]){
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
    for(int ic = 0; ic < ncells; ++ic){
        
        vector<int> fids = idFacesFromCell[ic];
        vector<int> nsigns = cnsigns[ic];
        
        for(int j = 0; j < fids.size(); ++j){
            int nsign = nsigns[j];
            int idf = fids[j];

            double fnx = fnormals[idf].first, fny = fnormals[idf].second;

            // b_mom_x = -Σ_f p_f * A_f * n_{x,f}
            b_mom_x[ic] = b_mom_x[ic] - p_face[idf] * fnx * nsign * flengths[idf];
            
            // b_mom_y = -Σ_f p_f * A_f * n_{y,f}
            b_mom_y[ic] = b_mom_y[ic] - p_face[idf] * fny * nsign * flengths[idf];
        }

        // + Não esquecer, ap já está divido por λ, por isso multiplica por (1-λ) ao invés de (1-λ)/λ.
        // & Adiciona contribuição de regularização: (1-λ) * aP * u_O^n 
        b_mom_x[ic] = b_mom_x[ic] + regularization * ap[ic] * uc[ic];
        // & Adiciona contribuição de regularização: (1-λ) * aP * v_O^n 
        b_mom_y[ic] = b_mom_y[ic] + regularization * ap[ic] * vc[ic];
    }

    A.setFromTriplets(triplets.begin(), triplets.end());
}

Eigen::VectorXd ReFumSolver::linear_upwind(Eigen::VectorXd b, vector<pair<BoundaryType, double>> boundary) {
    int ncells = mesh->get_ncells();

    for (int ic = 0; ic < ncells; ++ic) {
        vector<int> fids = idFacesFromCell[ic];
        vector<int> nsigns = cnsigns[ic];
        pair<double,double> centroidP = ccentroids[ic];
        
        for (int j = 0; j < fids.size(); ++j) {
            int idf = fids[j];
            int nsign = nsigns[j]; 
            
            // centro & normal corrigida
            pair<double,double> middleFace = fmiddles[idf];
            pair<double,double> normal = fnormals[idf];
            pair<double,double> normal_corrected = make_pair(normal.first * nsign, normal.second * nsign);

            double G_f = mdotf[idf] * nsign;

            if (fboundaryfaces[idf]) {
                
                if (boundary[idf].first == DIRICHLET) {
                    if (G_f > 0) {
                        // upwind = célula interna (P)
                        pair<double,double> deltaR = {middleFace.first - centroidP.first,
                                                      middleFace.second - centroidP.second};
                        double GdDotDeltaR = gradients(ic,0) * deltaR.first
                                           + gradients(ic,1) * deltaR.second;
                        
                        b[ic] -= G_f*GdDotDeltaR;
                    } 
                } 
                else if (boundary[idf].second == NEUMANN) {
                    continue;
                }
            } 
            else {
                // face interna
                if (G_f > 0) {                    
                    pair<double,double> deltaR = {middleFace.first - centroidP.first,
                                                  middleFace.second - centroidP.second};
                    double GdDotDeltaR = gradients(ic,0) * deltaR.first
                                       + gradients(ic,1) * deltaR.second;
                    
                    b[ic] -= G_f*GdDotDeltaR;
                } else {
                    int nb = get_neighbor(flftcs[idf], ic);
                    pair<double,double> centroidN = ccentroids[nb];
                    pair<double,double> deltaR = {middleFace.first - centroidN.first,
                                                  middleFace.second - centroidN.second};
                    double GdDotDeltaR = gradients(nb,0) * deltaR.first
                                       + gradients(nb,1) * deltaR.second;
                    
                    b[ic] -= G_f*GdDotDeltaR;
                }
            }
        }   
    }

    return b;
}

void ReFumSolver::solve_x_mom(){
    // Solve A x = b
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::DiagonalPreconditioner<double>> solver;
    solver.compute(A);
    solver.setTolerance(1e-4);

    Eigen::VectorXd aux = b_mom_x;
    for(int i = 0; i < 10; ++i){
        uc = solver.solve(aux); // resolve com valor atual de aux
        reconstruct_gradients(uc, u_face, u_boundary); // reconstroi e atualiza gradients.
        aux = cross_diffusion(b_mom_x, u_boundary); // att aux com cd.
        aux = linear_upwind(aux, u_boundary);
    }
    uc = solver.solve(aux);
}

void ReFumSolver::solve_y_mom(){
    // Solve A x = b 
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::DiagonalPreconditioner<double>> solver;
    solver.compute(A);
    solver.setTolerance(1e-4);

    Eigen::VectorXd aux = b_mom_y;
    for(int i = 0; i < 10; ++i){
        vc = solver.solve(aux); // resolve com 
        reconstruct_gradients(vc, v_face, v_boundary);
        aux = cross_diffusion(b_mom_y, v_boundary); // att aux com cd
        aux = linear_upwind(aux, v_boundary);
    }
    vc = solver.solve(aux);
}

/**
 * * Realiza a interpolação de momento para obter a velocidade nas faces. (PWIM)
 */
void ReFumSolver::face_velocity(){
    int nfaces = mesh->get_nedges();
    
    for(int idf = 0; idf < nfaces; ++idf){
        
        pair<int,int> nodes_share_face = flftcs[idf];
        int c1 = nodes_share_face.first;
        int c2 = nodes_share_face.second;

        if(fboundaryfaces[idf]){
            continue;
            /* 
            if(u_boundary[idf].first == DIRICHLET) 
                 continue; // Valor já está definido (wall ou inlet) e não é necessário alterá-lo.
             else{
                 As velocidades u_face e v_face quando for neumann, irei pegar do vizinho lá no momento, então aqui não precisa alterar nada.
                 continue;
            }*/
        }
        else{
            // * Face interior.
            double fnx = fnormals[idf].first, fny = fnormals[idf].second;

            // interpolação linear
            double velf_x = wf[idf]*uc[c1] + (1-wf[idf]) * uc[c2];
            double velf_y = wf[idf]*vc[c1] + (1-wf[idf]) * vc[c2];
            double vel_f_i = velf_x * fnx + velf_y * fny;

            double v0dp0_x = 0;
            double v0dp0_y = 0;

            vector<int> facesOfc1 = idFacesFromCell[c1];
            vector<int> nsignsc1 = cnsigns[c1];
            for(int j = 0; j < facesOfc1.size(); j++){
                // faces of cell O
                int faceC1 = facesOfc1[j];
                int nsign = nsignsc1[j];
                double fc1nx = fnormals[faceC1].first, fc1ny = fnormals[faceC1].second;
                
                // Σ p_f * A_f * n_x,f onde f são as faces da celula 1
                v0dp0_x = v0dp0_x + p_face[faceC1] * flengths[faceC1] * fc1nx * nsign;
                // Σ p_f * A_f * n_y,f onde f sao as faces da celula 1
                v0dp0_y = v0dp0_y + p_face[faceC1] * flengths[faceC1] * fc1ny * nsign;
            }

            double v1dp1_x = 0;
            double v1dp1_y = 0;

            vector<int> facesOfc2 = idFacesFromCell[c2];
            vector<int> nsignsc2 = cnsigns[c2];
            for(int j = 0; j < facesOfc2.size(); j++){
                // faces of cell N, contas similares ao caso de cima
                int faceC2 = facesOfc2[j];
                int nsign = nsignsc2[j];
                double fc2nx = fnormals[faceC2].first, fc2ny = fnormals[faceC2].second;

                v1dp1_x = v1dp1_x + p_face[faceC2] * flengths[faceC2] * fc2nx * nsign;
                v1dp1_y = v1dp1_y + p_face[faceC2] * flengths[faceC2] * fc2ny * nsign;
            }
            
            // continua contas do PWIM
            velf_x = wf[idf] * v0dp0_x/ap[c1] + (1-wf[idf]) * v1dp1_x/ap[c2];
            velf_y = wf[idf] * v0dp0_y/ap[c1] + (1-wf[idf]) * v1dp1_y/ap[c2];
            
            double velf_p = velf_x * fnx + velf_y * fny;
            
            double vdotn = vel_f_i + velf_p - (wf[idf]*careas[c1]/ap[c1] + (1-wf[idf])*careas[c2]/ap[c2]) * (pc[c2] - pc[c1])/fdfs[idf];
            
            // no fim, tem-se o fluxo de massa da face interpolado.
            mdotf[idf] = vdotn * rho * flengths[idf];
        }
    }
}

/**
 * * Calcula a correção da pressão (p')
 */
void ReFumSolver::solve_pp(bool sing_matrix){
    // *Calcula termo fonte: -Σ m_f
    int ncells = mesh->get_ncells();
    for(int ic = 0; ic < ncells; ++ic){
            
        // * warm-up
        vector<int> fids = idFacesFromCell[ic];
        vector<int> nsigns = cnsigns[ic];
        b_pc[ic] = 0; // zera para nn ter conteudo da anterior.
        
        for(int j = 0; j < fids.size(); ++j){
            //* warm-up
            int idf = fids[j];
            int nsign = nsigns[j];
            
            if(u_boundary[idf].first == NEUMANN){
                // Análogo ao caso do momento, é preciso usar as velocidades das celulas vizinhas para encontrar o mf.
                pair<double,double> normal_corrected = {
                        fnormals[idf].first * nsign,
                        fnormals[idf].second * nsign
                };
                double mf_outlet = rho * (uc[ic] * normal_corrected.first + vc[ic] * normal_corrected.second) * flengths[idf];
                b_pc[ic] = b_pc[ic] - mf_outlet;
            }else{
                // Quando for inlet ou wall (dirichlet) => usa o valor que já tem armazenado.
                b_pc[ic] = b_pc[ic] - mdotf[idf] * nsign;
            }
        }
    }
        
   // Calcula coeficientes da A
   triplets.clear();
   diags.setZero();
   for(int ic = 0; ic < ncells; ++ic){

        vector<int> fids = idFacesFromCell[ic];
        vector<int> nsigns = cnsigns[ic];

        for(int j = 0; j < fids.size(); ++j){ // loop over faces of cell
            int sign = nsigns[j];
            int idf = fids[j];

            if(fboundaryfaces[idf]){
                if(u_boundary[idf].first == DIRICHLET){
                    // se for wall, a correção do fluxo de massa é zero, então não tem nada para fazer.
                    // se for inlet, a correção do fluxo de massa é zero também, visto que já foi especificado, então não há nada a fazer.
                    continue;
                }else{
                    // se for neumann (outlet), então...
                    // correção do fluxo de massa em geral = rho_f * Af * (V_P/aP + V_N/aN) * (p'_O - p'_N)/delta_f
                    // tomar-se então: rho_f * Af * (V_P/aP)(p'_O - p'_b)/delta_f. Mas, p'b = 0 pois a pressão é dada.
                    // => rho_f * Af * (V_P/aP) * (p'_O)/delta_f 
                    // A_pc(ic, ic) = A_pc(ic,ic) + (rho * flengths[idf] * (careas[ic]/ap[ic])) / fdfs[idf];
                    diags[ic] = diags[ic] + (rho * flengths[idf] * (careas[ic]/ap[ic])) / fdfs[idf];
                }
            }
            else{

                // faces interiores
                int N = get_neighbor(flftcs[idf], ic);
                
                // A_pc(ic,ic) = A_pc(ic,ic) + rho * flengths[idf] * (wf[idf] * careas[ic]/ap[ic] + (1-wf[idf]) * careas[N]/ap[N])/fdfs[idf];
                diags[ic] = diags[ic] + rho * flengths[idf] * (wf[idf] * careas[ic]/ap[ic] + (1-wf[idf]) * careas[N]/ap[N])/fdfs[idf];

                // A_pc(ic, N) = -rho * flengths[idf] * (wf[idf] * careas[ic]/ap[ic] + (1-wf[idf]) * careas[N]/ap[N])/fdfs[idf];
                triplets.push_back(Eigen::Triplet<double>(ic, N, -rho * flengths[idf] * (wf[idf] * careas[ic]/ap[ic] + (1-wf[idf]) * careas[N]/ap[N])/fdfs[idf] ));
            }
        }
        triplets.push_back(Eigen::Triplet<double>(ic, ic, diags[ic]));
    }

    A.setFromTriplets(triplets.begin(), triplets.end());
    if(sing_matrix){
        int row = 0;
        for(int j = 0; j < A.outerSize(); ++j) {
            for(Eigen::SparseMatrix<double>::InnerIterator it(A,j); it; ++it) {
                if(it.row() == row && it.col() != row) it.valueRef() = 0.0;
            }
        }

        A.coeffRef(row, row) = 1.0;
        b_pc[row] = 0.0;
    }
    
    // Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::DiagonalPreconditioner<double>> solver;
    solver.setMaxIterations(100);
    solver.setTolerance(1e-3);
    solver.compute(A);
    pcorr = solver.solve(b_pc);
}

/*
* * Corrige velocidades.
*/
void ReFumSolver::uv_correct() {
    int nfaces = mesh->get_nedges();
    int ncells = mesh->get_ncells();

    // Interpola pfcorr nas faces
    for(int idf = 0; idf < nfaces; idf++){
        
        pair<int,int> id_nodes_share_face = flftcs[idf];
        int ic1 = id_nodes_share_face.first;
        int ic2 = id_nodes_share_face.second;

        if(fboundaryfaces[idf]){
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
    
    // + convergência (!!!) 
    uc_aux = uc; // * salva valor antigo
    vc_aux = vc; // * salva valor antigo
    
    for(int ic = 0; ic < ncells; ic++){
        
        ucorr[ic] = 0;
        vcorr[ic] = 0;
        
        vector<int> fids = idFacesFromCell[ic];
        vector<int> nsigns = cnsigns[ic];
        for(int j = 0; j < fids.size(); j++){
            int idf = fids[j];
            int nsign = nsigns[j];

            ucorr[ic] = ucorr[ic] + pfcorr[idf] * fnormals[idf].first * nsign * flengths[idf];
            vcorr[ic] = vcorr[ic] + pfcorr[idf] * fnormals[idf].second * nsign * flengths[idf];
        }
        ucorr[ic] = -ucorr[ic]/ap[ic];
        vcorr[ic] = -vcorr[ic]/ap[ic];
        
        uc[ic] = uc[ic] + ucorr[ic];
        vc[ic] = vc[ic] + vcorr[ic];
    }

    // Corrigindo valores das velocidades nas faces
    for(int idf = 0; idf < nfaces; idf++){
        
        if(fboundaryfaces[idf]){
            if(u_boundary[idf].first == DIRICHLET) 
                continue; // nada a fazer.
            else{
                // nada a fazer, pois não será utilizado, pois sempre pego o vizinho da célula!
            }
        }
        else{
            // interior.
            pair<int,int> nodes_share_face = flftcs[idf];
            int c1 = nodes_share_face.first;
            int c2 = nodes_share_face.second;

            double coeff = wf[idf]*careas[c1]/ap[c1] + (1-wf[idf]) * careas[c2]/ap[c2];
            mdotfcorr[idf] = rho*coeff*flengths[idf]*(pcorr[c1]-pcorr[c2])/fdfs[idf];
            mdotf[idf] = mdotf[idf] + mdotfcorr[idf];
        }
    }

};

void ReFumSolver::pres_correct(double lambda_p) {
    int ncells = mesh->get_ncells();
    
    pc_aux = pc; // * salva valores antigos
    for(int ic = 0; ic < ncells; ic++){
        pc[ic] = pc[ic] + lambda_p * pcorr[ic];
    }
}

// *=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=


void ReFumSolver::STEADY_SIMPLE(string problem, string filepath, int num_simple_iterations, double lambda_uv, double lambda_p, bool pressure_correction_flag){
    this->compute_bcs_repeat(); // aplica BC nos vetores
    double err_u = 1;
    double err_v = 1;
    double err_p = 1;

    // itera pelo número de iterações passado pelo usuário
    while(err_u > 1e-6 || err_v > 1e-6 || err_p > 1e-6){
        
        // cout << "# Calculando A_mom, b_mom_x e b_mom_y\n";
        mom_links_and_sources(lambda_uv);
        
        // cout << "Resolvendo para encontrar uc\n";
        solve_x_mom();
        
        // cout << "Resolvendo para encontrar uv\n";
        solve_y_mom();
        
        // cout << "# Calculando velocidade nas faces {u_f e v_f}\n";
        face_velocity();
        
        // cout << "# Calculando correção na pressão (p')\n";
        solve_pp(pressure_correction_flag); // lid chama com true, backward facing step chama com false
        
        // cout << "# Atualiza velocidades...\n";
        uv_correct();
        
        // cout << "# Atualiza pressão....\n";
        pres_correct(lambda_p);

        err_u = (uc - uc_aux).lpNorm<Eigen::Infinity>();
        err_v = (vc - vc_aux).lpNorm<Eigen::Infinity>();
        err_p = (pc - pc_aux).lpNorm<Eigen::Infinity>();
        cout << "\n||du||_inf = " << err_u << " ||dv||_inf = " << err_v << " ||dp||_inf = " << err_p << endl;

    }
    
    this->export_vtk(filepath + problem + ".vtk");

    this->calculate_exact_solution_and_compare();
}

void ReFumSolver::TRANSIENTE_SIMPLE(string problem, string filepath, int num_simple_iterations, double lambda_uv, double lambda_p, int n_steps, double tf, bool pressure_correction_flag){
    
    // calculando o dt (ASSUMINDO Q condição inicial é t = 0.)
    double t = 0;
    this->dt = (tf - 0)/n_steps;
    
    int cont = 0; // ajuda a salvar os arquivos.

    while(t <= tf){    
        
        // zera todos para a próxima iteração de tempo
        uc.setZero();
        vc.setZero();
        pc.setZero();
        u_face.setZero();
        v_face.setZero();
        p_face.setZero();
        mdotf.setZero();
        ap.setZero();
        mdotfcorr.setZero();
        pfcorr.setZero();
        pcorr.setZero();
        ucorr.setZero();
        vcorr.setZero();
        // -----------------------------------------------

        this->compute_bcs_repeat(); // aplica BC nos vetores
        
        for(int j = 0; j < num_simple_iterations; j++){
            
            // cout << "# Calculando A_mom, b_mom_x e b_mom_y\n";
            mom_links_and_sources(lambda_uv);
            
            // cout << "Resolvendo para encontrar uc\n";
            solve_x_mom();
            
            // cout << "Resolvendo para encontrar uv\n";
            solve_y_mom();
            
            // cout << "# Calculando velocidade nas faces {u_f e v_f}\n";
            face_velocity();
            
            // cout << "# Calculando correção na pressão (p')\n";
            solve_pp(true); // lid chama com true, backward facing step chama com false
            
            // cout << "# Atualiza velocidades...\n";
            uv_correct();
            
            // cout << "# Atualiza pressão....\n";
            pres_correct(lambda_p);
        }

        // faz variaveis na proxima iteração começarem como as antigas u(n) = u(n-1)
        uc_old = uc;
        vc_old = vc;
        pc_old = pc;
        
        string fn = filepath + problem + to_string(cont) + ".vtk";
        export_vtk(fn);

        t = t + dt; // atualiza o tempo
        cont = cont + 1; // atualiza contador de snapshots
        
        progress_bar(cont, 50, n_steps);
    }
    cout << endl;
}

void ReFumSolver::calculate_exact_solution_and_compare(){
    int ncells = mesh->get_ncells();
    Eigen::VectorXd u_exact(ncells);
    u_exact.setZero();
    Eigen::VectorXd v_exact(ncells);
    v_exact.setZero();
    Eigen::VectorXd p_exact(ncells);
    p_exact.setZero();

    vector<Cell*> cells = mesh->get_cells();
    double lambda = -1.81009812;
    for(int i = 0; i < ncells; ++i){
        Cell* c = cells[i];
        pair<double,double> centroid = c->get_centroid();
        double xc = centroid.first; double yc = centroid.second;
        u_exact[i] = 1 - exp(lambda*xc) * cos(2*M_PI*yc);
        v_exact[i] = (lambda/(2*M_PI))*exp(lambda*xc)*sin(2*M_PI*yc);
        p_exact[i] = 0.5 * (1 - exp(2 * lambda * xc));
    }
    cout << "==*===*===*==*===*===*==*===*===*==*===*===*==*===*===*\n";
    cout << "\nDiferença entre u_exact e u:" << (u_exact - uc).lpNorm<Eigen::Infinity>() << endl;
    cout << "Diferença entre v_exact e v:" << (v_exact - vc).lpNorm<Eigen::Infinity>() << endl;
    cout << "Diferença entre p_exact e p:" << (p_exact - pc).lpNorm<Eigen::Infinity>() << endl;
}


// *=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=*=**=*=


void ReFumSolver::export_vtk(string filename) {

    ofstream vtk_file(filename);

    if (!vtk_file.is_open()) {
        std::cerr << "ERROR: failed to create file '" << filename << "'.\n";
        exit(1);
    }

    // ---------------- Informações padrão ----------------
    vtk_file << "# vtk DataFile Version 3.0\n";
    vtk_file << "Malha 2D com triângulos: velocidade e pressão\n";
    vtk_file << "ASCII\n";
    vtk_file << "DATASET UNSTRUCTURED_GRID\n\n";

    // ---------------- Pontos ----------------
    vector<Node*> nodes = this->mesh->get_nodes();
    vtk_file << "POINTS " << nodes.size() << " double" << endl;
    for (auto &node : nodes) {
        vtk_file << node->x << " " << node->y << " 0\n";
    }
    vtk_file << endl;

    // ---------------- Células ----------------
    vector<Cell*> cells = mesh->get_cells();
    vtk_file << "CELLS " << cells.size() << " ";
    int sum = 0;
    for (auto &cell : cells) {
        sum += cell->get_nodes().size() + 1; // +1 para o número de nós no VTK
    }
    vtk_file << sum << endl;

    for (auto &cell : cells) {
        vector<Node*>& nodesOfCell = cell->get_nodes();
        vtk_file << nodesOfCell.size();
        for (auto &n : nodesOfCell)
            vtk_file << " " << n->id;
        vtk_file << endl;
    }
    vtk_file << endl;

    // ---------------- Tipos de células ----------------
    vtk_file << "CELL_TYPES " << cells.size() << endl;
    for (auto &cell : cells) {
        if (cell->cellType == 2)
            vtk_file << "5\n"; // triângulo
        else if (cell->cellType == 3)
            vtk_file << "9\n"; // quadrilátero
    }
    vtk_file << endl;

    // ---------------- Dados das células ----------------
    vtk_file << "CELL_DATA " << cells.size() << endl;

    // Velocidade
    vtk_file << "VECTORS velocity double\n";
    for (size_t i = 0; i < cells.size(); ++i) {
        vtk_file << uc[i] << " " << vc[i] << " 0\n";
    }
    vtk_file << endl;

    // Pressão
    vtk_file << "SCALARS pressure double 1\n";
    vtk_file << "LOOKUP_TABLE default\n";
    for (size_t i = 0; i < cells.size(); ++i) {
        vtk_file << pc[i] << endl;
    }

    vtk_file.close();
}
