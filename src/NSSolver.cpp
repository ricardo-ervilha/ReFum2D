#include "NSSolver.h"
#include "Mesh.h"
#include "Edge.h"
#include "Cell.h"
#include "BoundaryCondition.h"

NSSolver::NSSolver(Mesh *mesh, float mu, float rho, float source_x, float source_y){
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

    // Utilizados na parte de calcular momento (u*, v*)
    this->A_mom = arma::mat(ncells, ncells);
    this->b_mom = arma::vec(ncells, arma::fill::zeros);

    // Utilizados para calcular a correção da pressão (p')
    this->A_pc = arma::mat(ncells, ncells);
    this->b_pc = arma::vec(ncells, arma::fill::zeros);
    
    // velocidades corrigidas nos centroides
    this->u = arma::vec(ncells, arma::fill::zeros);
    this->v = arma::vec(ncells, arma::fill::zeros);

    // velocidades aproximadas nos centroides
    this->u_star = arma::vec(ncells, arma::fill::zeros);
    this->v_star = arma::vec(ncells, arma::fill::zeros);

    // velocidades nas faces interpoladas pela interpolação de momento
    this->u_face = arma::vec(nfaces, arma::fill::zeros);
    this->v_face = arma::vec(nfaces, arma::fill::zeros);

    this->mdotf = arma::vec(nfaces, arma::fill::zeros);

    this->u_face_prime = arma::vec(nfaces, arma::fill::zeros);
    this->v_face_prime = arma::vec(nfaces, arma::fill::zeros);

    // variáveis relacionadas a pressão
    this->p = arma::vec(ncells, arma::fill::zeros);
    this->p_prime = arma::vec(ncells, arma::fill::zeros);
    this->p_star = arma::vec(ncells, arma::fill::zeros);

    // coeficientes de Difusão e Convecção
    this->Df = arma::vec(nfaces, arma::fill::zeros);
    this->Gf = arma::vec(nfaces, arma::fill::zeros);

    this->aP = arma::vec(ncells, arma::fill::zeros);

    // =======================================================================================

    // TODO: Preenche boundary top do lid-driven cavity flow.
    vector<Edge*> faces = mesh->get_edges();
    for(int i = 0; i < nfaces; i++){
        Edge* face = faces[i];
        if(face->from->y == 0.01 && face->to->y == 0.01){
            // Top 
            u_face[face->id] = 0.01; 
        }
    }
}   

arma::vec jacobi(const arma::mat& A, const arma::vec& b, int max_iter = 100000, double tol = 1e-6) {
    int n = A.n_rows;
    arma::vec x = arma::zeros<arma::vec>(n);     // chute inicial
    arma::vec x_new = x;

    for (int k = 0; k < max_iter; k++) {
        for (int i = 0; i < n; i++) {
            double sigma = 0.0;
            for (int j = 0; j < n; j++) {
                if (j != i)
                    sigma += A(i, j) * x(j);
            }
            x_new(i) = (b(i) - sigma) / A(i, i);
        }

        // critério de parada
        if (norm(x_new - x, "inf") < tol)
            break;

        if(k % 1000 == 0)
            cout << norm(x_new - x, "inf") << endl;

        x = x_new;
    }

    return x_new;
}

NSSolver::~NSSolver(){
    // nada.
}

/**
 * * Calcula a matriz A da equação de momentum. Igual para u e v!
 */
void NSSolver::calculate_A_mom(){
    vector<Cell*> cells = mesh->get_cells();
    
    for(int i = 0; i < cells.size(); ++i){
        // * warm-up
        Cell* c = cells[i];
        vector<int> &nsigns = c->get_nsigns();
        vector<Edge*> faces = c->get_edges();
        
        double a_P = 0; // * sum_f (max(0, ṁ_f) + D_f)
       
        for(int j = 0; j < faces.size(); ++j){
            // * warm-up
            Edge* face = faces[j];
            int nsign = nsigns[j];
            pair<double, double> &normal = face->get_normal();
            pair<double, double> normal_corrected = make_pair(normal.first * nsign, normal.second * nsign);
            
            double Df = (mu * face->get_length()) / face->get_df();
            double a_N = min(0.0, mdotf[face->id] * nsign) - Df;
            
            
            if(!face->is_boundary_face()){
                // * contabiliza a_N na matriz A (off-diagonal)
                int N = get_neighbor(face->get_link_face_to_cell(), c->id);
                this->A_mom(c->id, N) += a_N;
                a_P += max(0.0, mdotf[face->id] * nsign) + Df;
            }else{
                // * poderá ser utilizado o valor da face de contorno para contabilizar como termo fonte.
                a_P += Df;
            }
        }

        // * Contabiliza o a_P na diagonal
        this->A_mom(c->id, c->id) += a_P;

        // * salva o a_P para uso posterior
        aP[c->id] = a_P; // sobrescreve
    }
}

/**
 * * Calcula o vetor b do sistema para x.
 */
void NSSolver::calculate_b_mom_x(){
    vector<Cell*> cells = mesh->get_cells();
    for(int i = 0; i < cells.size(); ++i){
        // * Warm-up
        Cell* c = cells[i];
        vector<Edge*> faces = c->get_edges();
        vector<int> &nsigns = c->get_nsigns();

        // * armazena: \sum_{b(f)} a_f u_f
        double boundary_contribution = 0;

        // * - sum_f (p_f) * n_xf * A_f
        double pressure_contribution = 0;
        
        for(int j = 0; j < faces.size(); ++j){
            int nsign = nsigns[j];
            Edge* face = faces[j];
            pair<double, double> &normal = face->get_normal();
            pair<double, double> normal_corrected = make_pair(normal.first * nsign, normal.second * nsign);

            if(face->is_boundary_face()){
                double Df = (mu * face->get_length()) / face->get_df();
               
                boundary_contribution += Df * u_face[face->id] - mdotf[face->id] * nsign * u_face[face->id]; // * a_F u_F

                // & Tratamento da pressão
                // estou usando série de taylor: phi_P + grad(Phi)_N * d_Pf;
                // ou seja, o valor de phi na face de contorno é o valor na celula que está perto dela + gradiente na face * distancia entre eles. Como no problema o gradiente é sempre zero, então o phi na face é o phi do vizinho.
                double pf = this->p_star[c->id];
                double nxf = normal_corrected.first;
                pressure_contribution += pf * nxf * face->get_length();
            }else{
                // * Nada a fazer em termos de contorno, pois já contabilizou na A.
                
                // & Tratamento da pressão
                int N = get_neighbor(face->get_link_face_to_cell(), c->id);
                // pressão na face obtida pela média entre os vizinhos
                double pf = 0.5 * (this->p_star[c->id] + this->p_star[N]);
                double nxf = normal_corrected.first;
                pressure_contribution += pf * nxf * face->get_length();
            }
        }
        
        // * Termo fonte
        // double source_contribution = source_x * c->get_area(); // TODO: AINDA NÃO USA

        this->b_mom[c->id] += -pressure_contribution + boundary_contribution;
    }
}

/**
 * * Calcula o vetor b do sistema para y.
 */
void NSSolver::calculate_b_mom_y(){
    vector<Cell*> cells = mesh->get_cells();
    for(int i = 0; i < cells.size(); ++i){
        // * Warm-up
        Cell* c = cells[i];
        vector<Edge*> faces = c->get_edges();
        vector<int> &nsigns = c->get_nsigns();

        // * armazena: \sum_{b(f)} a_f v_f
        double boundary_contribution = 0;

        // * - sum_f (p_f) * n_yf * A_f
        double pressure_contribution = 0;
        
        for(int j = 0; j < faces.size(); ++j){
            int nsign = nsigns[j];
            Edge* face = faces[j];
            pair<double, double> &normal = face->get_normal();
            pair<double, double> normal_corrected = make_pair(normal.first * nsign, normal.second * nsign);

            if(face->is_boundary_face()){
                double Df = (mu * face->get_length()) / face->get_df();
               
                boundary_contribution += Df * v_face[face->id] - mdotf[face->id] * nsign * v_face[face->id]; // * a_F v_F

                // & Tratamento da pressão
                // estou usando série de taylor: phi_P + grad(Phi)_N * d_Pf;
                // ou seja, o valor de phi na face de contorno é o valor na celula que está perto dela + gradiente na face * distancia entre eles. Como no problema o gradiente é sempre zero, então o phi na face é o phi do vizinho.
                double pf = this->p_star[c->id];
                double nyf = normal_corrected.second;
                pressure_contribution += pf * nyf * face->get_length();
            }else{
                // * Nada a fazer em termos de contorno, pois já contabilizou na A.
                
                // & Tratamento da pressão
                int N = get_neighbor(face->get_link_face_to_cell(), c->id);
                // pressão na face obtida pela média entre os vizinhos
                double pf = 0.5 * (this->p_star[c->id] + this->p_star[N]);
                double nyf = normal_corrected.second;
                pressure_contribution += pf * nyf * face->get_length();
            }
        }
        
        // * Termo fonte
        // double source_contribution = source_y * c->get_area(); // TODO: AINDA NÃO USA

        this->b_mom[c->id] += -pressure_contribution + boundary_contribution;
    }
}

/**
 * * Resolve as equações de momento para obter u* e v*!
 */
void NSSolver::calculate_momentum(){
    cout << "Calculando as equações de momento.\n";

    // * Zera para desconsiderar contas anteriores feitas.
    this->A_mom.zeros();
    this->b_mom.zeros();

    // * Preenche A
    cout << "Calculando os coeficientes da matriz A.\n";
    calculate_A_mom();

    // * Preenche b para o momento em x
    cout << "Calculando os coeficientes do vetor b de x-mom.\n";
    calculate_b_mom_x();

    // cout << A_mom << endl;
    // cout << b_mom << endl;

    cout << "Calculando u_star.\n";
    this->u_star = jacobi(this->A_mom, this->b_mom); // + Resolve

    // * Zera pra calcular o do momento em y
    this->b_mom.zeros();

    cout << "Calculando os coeficientes do vetor b de y-mom.\n";
    calculate_b_mom_y();
    
    cout << "Calculando v_star.\n";
    this->v_star = jacobi(this->A_mom, this->b_mom); // + Resolve

    cout << "=============================================\n";
}

/**
 * * Realiza a interpolação de momento para obter a velocidade nas faces.
 */
void NSSolver::interpolate_momentum(){
    cout << "Interpolando a velocidade nas faces.\n";
    vector<Edge*> faces = mesh->get_edges();
    vector<Cell*> cells = mesh->get_cells();
    
    for(int i = 0; i < faces.size(); ++i){
        Edge* face = faces[i];
        if(!face->is_boundary_face()){ // * Não pode ser boundary face.
            pair<int,int> id_cells_share_face = face->get_link_face_to_cell();
            int P = id_cells_share_face.first;
            int N = id_cells_share_face.second;
            
            // * u_f = (u_P + u_N)/2 => média
            double u_bar = 0.5 * (u_star[P] + u_star[N]);
            double v_bar = 0.5 * (v_star[P] + v_star[N]);
            double U_bar_dot_normal = u_bar * face->get_normal().first + v_bar * face->get_normal().second;

            double V0dp0_x = 0;
            double V0dp0_y = 0;
            vector<int> nsignsP = cells[P]->get_nsigns();
            vector<Edge*> facesOfP = cells[P]->get_edges();
            for(int j = 0; j < facesOfP.size(); ++j){
                pair<double, double> normal = facesOfP[j]->get_normal();
                pair<double, double> normal_corrected = make_pair(normal.first * nsignsP[j], normal.second * nsignsP[j]);
                double pf;
                if(facesOfP[j]->is_boundary_face()){
                    pf = this->p_star[P];
                }else{
                    int Neighbor = get_neighbor(facesOfP[j]->get_link_face_to_cell(), P);
                    pf = 0.5 * (this->p_star[P] + this->p_star[Neighbor]);
                }
                V0dp0_x += pf * facesOfP[j]->get_length() * normal_corrected.first;
                V0dp0_y += pf * facesOfP[j]->get_length() * normal_corrected.second;
            }

            double V1dp1_x = 0;
            double V1dp1_y = 0;
            vector<int> nsignsN = cells[N]->get_nsigns();
            vector<Edge*> facesOfN = cells[N]->get_edges();
            for(int j = 0; j < facesOfN.size(); ++j){
                pair<double, double> normal = facesOfN[j]->get_normal();
                pair<double, double> normal_corrected = make_pair(normal.first * nsignsN[j], normal.second * nsignsN[j]);
                double pf;
                if(facesOfN[j]->is_boundary_face()){
                    pf = this->p_star[N];
                }else{
                    int Neighbor = get_neighbor(facesOfN[j]->get_link_face_to_cell(), N);
                    pf = 0.5 * (this->p_star[N] + this->p_star[Neighbor]);
                }
                V1dp1_x += pf * facesOfN[j]->get_length() * normal_corrected.first;
                V1dp1_y += pf * facesOfN[j]->get_length() * normal_corrected.second;
            }

            double velf_x = 0.5 * V0dp0_x/aP[P] + 0.5 * V1dp1_x /aP[N];
            double velf_y = 0.5 * V0dp0_y/aP[P] + 0.5 * V1dp1_y /aP[N];
            double velf_p = velf_x * face->get_normal().first +  velf_y * face->get_normal().second;

            double term_test = velf_p - 0.5 * (cells[P]->get_area()/aP[P] + cells[N]->get_area()/aP[N]) * ((p_star[N] - p_star[P])/face->get_df());

            double vdotn = U_bar_dot_normal + term_test;

            mdotf[face->id] = vdotn * rho * face->get_length(); // mass flow rate.

        }else{
            // * já tem valor bem definido.
        }
    }
    // for(int i = 0; i < mdotf.size(); i++){
    //     cout << "x: " << faces[i]->get_middle().first << "\ty: " << faces[i]->get_middle().second << "\tuf: " << mdotf[i] << endl;
    // }
    cout << "=============================================\n";
}

/**
 * * Calcula a correção da pressão (p')
 */
void NSSolver::pressure_correction_poisson(){
    cout << "Resolvendo a poisson para encontrar a correção da pressão.\n";
    
    // * zera pra próxima iteração
    A_pc.zeros();
    b_pc.zeros();
    p_prime.zeros();

    vector<Cell*> cells = mesh->get_cells();
    for(int i = 0; i < cells.size(); ++i){
        
        // * warm-up
        Cell* c = cells[i];
        vector<Edge*> faces = c->get_edges();
        vector<int> nsigns = c->get_nsigns();
        
        double a_P = 0; // * coeficiente do p'_C
        double b = 0; // * valor b do lado direito
        
        for(int j = 0; j < faces.size(); ++j){
            //* warm-up
            Edge* face = faces[j];
            int nsign = nsigns[j];

            if(face->is_boundary_face()){
                // * skip
            }else{
                // * contribue no vizinho da diagonal
                int N = get_neighbor(face->get_link_face_to_cell(), c->id);
                
                double dP = cells[c->id]->get_area()/aP[c->id];
                double dN = cells[N]->get_area()/aP[N];
                double a_N = 0.5 * (dP + dN) * (rho * face->get_length() / face->get_df());

                this->A_pc(c->id, N) += -a_N; // off-diagonal

                a_P += a_N; // acumula os ans
            }

            b += mdotf[face->id] * nsign; // acumula o fluxo de massa
        }
        
        // * diagonal
        A_pc(c->id, c->id) += a_P;
        b_pc(c->id) += -b; // * fonte
    }
    
    // + Resolvendo o problema de não ter valor prescrito de pressão...
    // A_pc.row(0).zeros(); // zera a linha
    // A_pc(0,0) = 1.0;
    // b_pc[0] = 0.0; // prescrevendo uma correção em um nó. Deixará de ser singular a matriz...
    
    this->p_prime = jacobi(A_pc, b_pc, 100);

    cout << "=============================================\n";
}



/**
 * * Corrige as variáveis...
 */
void NSSolver::correct_variables(double alpha_uv, double alpha_p){
    cout << "Corrigindo as variáveis.\n";
    
    // corrige velocidade nos centroides
    vector<Cell*> cells = mesh->get_cells();
    for(int i = 0; i < cells.size(); i++){
        Cell* c = cells[i];
        vector<int> nsigns = c->get_nsigns();
        vector<Edge*> faces = c->get_edges();

        double correction_u = 0;
        double correction_v = 0;
        for(int j = 0; j < faces.size(); j++){
            Edge* face = faces[j];
            int nsign = nsigns[j];
            pair<double, double> normal = face->get_normal();
            pair<double, double> normal_corrected = make_pair(normal.first * nsign, normal.second * nsign);

            double p_prime_f;
            if(face->is_boundary_face()){
                p_prime_f = this->p_prime[c->id];
            }else{
                int Neighbor = get_neighbor(face->get_link_face_to_cell(), c->id);
                p_prime_f = 0.5 * (this->p_prime[c->id] + this->p_prime[Neighbor]);
            }

            correction_u += p_prime_f * face->get_length() * normal_corrected.first;
            correction_v += p_prime_f * face->get_length() * normal_corrected.second;
        }

        u[c->id] = u_star[c->id] - alpha_uv * correction_u/aP[c->id];
        v[c->id] = v_star[c->id] - alpha_uv * correction_v/aP[c->id];
    }

    // corrige mdotf
    vector<Edge*> faces = mesh->get_edges();
    for(int i = 0; i < faces.size(); i++){
        Edge* face = faces[i];
        if(!face->is_boundary_face()){ // * SÓ ATUALIZA SE NÃO É FACE DE CONTORNO
            pair<int,int> id_cells_share_face = face->get_link_face_to_cell();
            int P = id_cells_share_face.first;
            int N = id_cells_share_face.second;
            double dP = cells[P]->get_area()/aP[P];
            double dN = cells[N]->get_area()/aP[N];
            double dF = 0.5 * (dP + dN);
    
            mdotf[face->id] = mdotf[face->id] - alpha_uv * rho * face->get_length() * dF * (p_prime[N] - p_prime[P]) / face->get_df();
        }
    }

    // atualiza pressão;
    for(int i = 0; i < cells.size(); i++){
        this->p[i] = this->p_star[i] + alpha_p * p_prime[i];
    }

    // update
    u_star = u;
    v_star = v;
    p_star = p;
}

void NSSolver::export_solution(string filename, char variable){
    arma::vec var;
    if(variable == 'u')
        var = this->u_star;
    else if(variable == 'v')
        var = this->v_star;
    else if(variable == 'p')
        var = this->p_star;

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
        vtk_file << var[i] << endl;
    }

    vtk_file.close();
}