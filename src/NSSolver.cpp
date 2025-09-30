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

    this->u_face_prime = arma::vec(nfaces, arma::fill::zeros);
    this->v_face_prime = arma::vec(nfaces, arma::fill::zeros);

    // variáveis relacionadas a pressão
    this->p = arma::vec(ncells, arma::fill::zeros);
    this->p_prime = arma::vec(ncells, arma::fill::zeros);
    this->p_star = arma::vec(ncells, arma::fill::zeros);

    // coeficientes de Difusão e Convecção
    this->Df = arma::vec(nfaces, arma::fill::zeros);
    this->Gf = arma::vec(nfaces, arma::fill::zeros);

    this->a_coeff = arma::vec(ncells, arma::fill::zeros);

    this->u_EXACT = arma::vec(ncells, arma::fill::zeros);
    this->v_EXACT = arma::vec(ncells, arma::fill::zeros);
    this->p_EXACT = arma::vec(ncells, arma::fill::zeros);

    // ============================== PREENCHENDO OS EXATOS =================================
    vector<Cell*> cells = mesh->get_cells();
    for(int i = 0; i < cells.size(); i++){
        u_EXACT[i] = cells[i]->get_centroid().first; // x
        v_EXACT[i] = - cells[i]->get_centroid().second; // - y
        p_EXACT[i] = - (pow(cells[i]->get_centroid().first, 2) + pow(cells[i]->get_centroid().second, 2))/2.0;
    }

    // =======================================================================================

    // preenchendo BOUNDARIES
    vector<Edge*> faces = mesh->get_edges();
    for(int i = 0; i < nfaces; i++){
        Edge* face = faces[i];
        if(face->from->y == 1 && face->to->y == 1){
            // Top 
            u_face[face->id] = 1; 
        }
        // if(face->from->y == 1 && face->to->y == 1){
        //     // Top 
        //     u_face[face->id] = face->get_middle().first; 
        //     v_face[face->id] = -1;
        // }
        // else if(face->from->y == 0 && face->to->y == 0){
        //     // Bottom wall
        //     u_face[face->id] = face->get_middle().first ;
        //     v_face[face->id] = 0.0;
        // }
        // else if(face->from->x == 0 && face->to->x == 0){
        //     // Left wall
        //     u_face[face->id] = 0.0;
        //     v_face[face->id] = -face->get_middle().second;
        // }
        // else if(face->from->x == 1 && face->to->x == 1){
        //     // Right wall
        //     u_face[face->id] = 1.0;
        //     v_face[face->id] = -face->get_middle().second; 
        // }
    }
}   

NSSolver::~NSSolver(){
    // nada.
}

void NSSolver::calculate_Df(){
    cout << "Calculando os valores de Df.\n";
    vector<Edge*> faces = this->mesh->get_edges();
    for(int i = 0; i < faces.size(); ++i){
        // * (μ_f*A_f)/d_CF
        Df[i] = mu * faces[i]->get_length() / faces[i]->get_df();
    }
    cout << "=============================================\n";
}

void NSSolver::calculate_Gf(){
    cout << "Calculando os valores de Gf.\n";
    vector<Edge*> faces = this->mesh->get_edges();
    for(int i = 0; i < faces.size(); ++i){
        // * ρ_f * A_f | nota: o restante que seria (v_f \cdot n_f) será calculado no momentum.
        Gf[i] = rho * faces[i]->get_length();
    }
    cout << "=============================================\n";
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
        
        double a_P = 0; // * sum_f (max(0, G_f) + D_f)
       
        for(int j = 0; j < faces.size(); ++j){
            // * warm-up
            Edge* face = faces[j];
            int nsign = nsigns[j];
            pair<double, double> &normal = face->get_normal();
            pair<double, double> normal_corrected = make_pair(normal.first * nsign, normal.second * nsign);
            
            // + Usando WMI
            // * v_f \cdot n_f 
            double UdotNormal = normal_corrected.first * u_face[face->id] + normal_corrected.second * v_face[face->id];

            double a_N = min(0.0, Gf[face->id] * UdotNormal) - Df[face->id];
            if(!face->is_boundary_face()){
                // * contabiliza a_N na matriz A (off-diagonal)
                int N = get_neighbor(face->get_link_face_to_cell(), c->id);
                this->A_mom(c->id, N) += a_N;
            }else{
                // * contabiliza a_N * u_F ou a_N * v_f no vetor b.
            }

            a_P += max(0.0, Gf[face->id] * UdotNormal) + Df[face->id];
        }

        // * Contabiliza o a_P na diagonal
        this->A_mom(c->id, c->id) += a_P;
        
        // * registra em a_coeff o a_P para ser utilizado depois em outros cálculos.
        a_coeff[c->id] = a_P;
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
        
        for(int j = 0; j < faces.size(); ++j){
            int nsign = nsigns[j];
            Edge* face = faces[j];
            pair<double, double> &normal = face->get_normal();
            pair<double, double> normal_corrected = make_pair(normal.first * nsign, normal.second * nsign);

            // + Usando o vetor de WMI
            double UdotNormal = normal_corrected.first * u_face[face->id] + normal_corrected.second * v_face[face->id];

            if(face->is_boundary_face()){
                double a_N = min(0.0, Gf[face->id] * UdotNormal) - Df[face->id];
                boundary_contribution += -a_N * u_face[face->id]; // * a_F u_F
            }else{
                // * Nada a fazer, pois já contabilizou na A.
            }
        }
        
        // * obtém (∇p*)_P
        pair<double, double> pressure_gradient = reconstruct_pressure_gradients(c, this->p_star);
        double pressure_contribution = - (pressure_gradient.first * c->get_area()); // pega \partial p / \partial x

        // * Termo fonte
        double source_contribution = source_x * c->get_area();

        this->b_mom[c->id] += pressure_contribution + source_contribution + boundary_contribution;
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
        
        for(int j = 0; j < faces.size(); ++j){
            int nsign = nsigns[j];
            Edge* face = faces[j];
            pair<double, double> &normal = face->get_normal();
            pair<double, double> normal_corrected = make_pair(normal.first * nsign, normal.second * nsign);

            // + Usando WMI
            double UdotNormal = normal_corrected.first * u_face[face->id] + normal_corrected.second * v_face[face->id];

            if(face->is_boundary_face()){
                double a_N = min(0.0, Gf[face->id] * UdotNormal) - Df[face->id];
                boundary_contribution += -a_N * v_face[face->id]; // * a_F * v_F
            }else{
                // * Nada a fazer, pois já contabilizou na A.
            }
        }
        
        // * obtém (∇p*)_P
        pair<double, double> pressure_gradient = reconstruct_pressure_gradients(c, this->p_star);
        double pressure_contribution = - (pressure_gradient.second * c->get_area()); // pega \partial p / \partial y

        // * Termo fonte
        double source_contribution = source_y * c->get_area();

        this->b_mom[c->id] += pressure_contribution + source_contribution + boundary_contribution;
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

    cout << "Calculando u_star.\n";
    this->u_star = arma::solve(this->A_mom, this->b_mom); // + Resolve

    // * Zera pra calcular o do momento em y
    this->b_mom.zeros();

    cout << "Calculando os coeficientes do vetor b de y-mom.\n";
    calculate_b_mom_y();
    
    cout << "Calculando v_star.\n";
    this->v_star = arma::solve(this->A_mom, this->b_mom); // + Resolve

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
            
            pair<int,int> adjacent_cells = face->get_link_face_to_cell();
            int cell_C = adjacent_cells.first;
            int cell_F = adjacent_cells.second;

            // * Interpolação da velocidade (v_bar) ---
            double u_avg = 0.5 * (u_star[cell_C] + u_star[cell_F]);
            double v_avg = 0.5 * (v_star[cell_C] + v_star[cell_F]);

            // * Vetor e_CF (normalizado) ---
            pair<double,double> dCF_vec = {
                cells[cell_F]->get_centroid().first - cells[cell_C]->get_centroid().first,
                cells[cell_F]->get_centroid().second - cells[cell_C]->get_centroid().second
            };
            double dCF = face->get_df();
            pair<double,double> eCF = { dCF_vec.first / dCF, dCF_vec.second / dCF }; // * normaliza

            // * Valor de d_f ---
            double dC = cells[cell_C]->get_area() / a_coeff[cell_C];
            double dF = cells[cell_F]->get_area() / a_coeff[cell_F];
            double d_f = 0.5 * (dC + dF);

            // * Termo (pF - pC)/dCF ---
            double pressure_diff = (p[cell_F] - p[cell_C]) / dCF;

            // gradiente de pressão interpolado na face
            auto gradP_C = reconstruct_pressure_gradients(cells[cell_C], this->p_star);
            auto gradP_F = reconstruct_pressure_gradients(cells[cell_F], this->p_star);
            pair<double,double> gradP_face = {
                0.5 * (gradP_C.first + gradP_F.first),
                0.5 * (gradP_C.second + gradP_F.second)
            };

            // * produto escalar gradPf ⋅ eCF
            double grad_dot_eCF = gradP_face.first * eCF.first + gradP_face.second * eCF.second;

            double correction = (pressure_diff - grad_dot_eCF);

            // * Correção final
            double u_corr = d_f * correction * eCF.first;
            double v_corr = d_f * correction * eCF.second;

            u_face[face->id] = u_avg - u_corr;
            v_face[face->id] = v_avg - v_corr;
        }
    }
    cout << "=============================================\n";
}

/**
 * * Função para reconstruir os gradientes quando necessário.
 */
pair<double,double> NSSolver::reconstruct_pressure_gradients(Cell *c, arma::vec phi){  

    vector<Cell*> cells = mesh->get_cells();

    vector<Edge*>& edgesOfCell = c->get_edges();
    pair<double,double>& centroidP = c->get_centroid();

    // criação da Matriz M e do vetor y.
    arma::mat M(edgesOfCell.size(), 2);
    arma::vec y(edgesOfCell.size());
        
    for(int j = 0; j < edgesOfCell.size(); j++){
        Edge* e = edgesOfCell[j];
        
        // # meio da face
        pair<double, double> &middleFace = e->get_middle();
        
        if(e->is_boundary_face()){
            // * considero distancia entre C e centro da face
            double dx = middleFace.first - centroidP.first;
            double dy = middleFace.second - centroidP.second;
            // TODO: Esse tratamento faz sentido ? Estou colocando zero pois não sei qual valor seria.
            M(j,0) = dx;
            M(j,1) = dy;
            y[j] = 0;
        }else{
            int nb = get_neighbor(e->get_link_face_to_cell(), c->id);
            pair<double,double>& centroidN = cells[nb]->get_centroid();

            double dx = centroidN.first - centroidP.first;
            double dy = centroidN.second - centroidP.second;
            
            // + [dx dy]
            M(j,0) = dx;
            M(j, 1) = dy;

            // + [u_N - u_P]
            y[j] = phi[nb] - phi[c->id];  
        }
    }
        
    // + Resolve sistema sobre-determinado M x = y.
    arma::vec x = arma::solve(M, y);
    double phix = x[0];
    double phiy = x[1];
    
    return make_pair(phix, phiy);
}

/**
 * * Calcula a correção da pressão (p')
 */
void NSSolver::pressure_correction_poisson(){
    cout << "Resolvendo a poisson para encontrar a correção da pressão.\n";
    
    // * zera pra próxima iteração
    A_pc.zeros();
    b_pc.zeros();

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
            pair<double,double> normal = face->get_normal();
            pair<double, double> normal_corrected = make_pair(normal.first * nsign, normal.second * nsign);

            double UdotNormal = normal_corrected.first * u_face[face->id] + normal_corrected.second * v_face[face->id];
            
            if(face->is_boundary_face()){
                // * contribue do lado direito como termo fonte
            }else{
                // * contribue na A
                int N = get_neighbor(face->get_link_face_to_cell(), c->id);
                
                double d_P = cells[c->id]->get_area()/a_coeff[c->id];
                double d_N = cells[N]->get_area()/a_coeff[N];
                double d_f = 0.5 * (d_P + d_N);

                // * valores de d_CF^x e y
                double dcfx = cells[N]->get_centroid().first - cells[c->id]->get_centroid().first;
                double dcfy = cells[N]->get_centroid().second - cells[c->id]->get_centroid().second;

                // * valores de S_f^x e y
                double sfx = face->get_length() * normal_corrected.first;
                double sfy = face->get_length() * normal_corrected.second;

                pair<double,double> dCF_vec = {
                    cells[N]->get_centroid().first - cells[c->id]->get_centroid().first,
                    cells[N]->get_centroid().second - cells[c->id]->get_centroid().second
                };
                double dCF = face->get_df();
                pair<double,double> eCF = { dCF_vec.first / dCF, dCF_vec.second / dCF }; // * normaliza

                double eCFdotNormal = eCF.first * normal_corrected.first + eCF.second * normal_corrected.second;

                double a_N = (d_f * eCFdotNormal * face->get_length())/face->get_df();
                
                A_pc(c->id, N) += -a_N;

                a_P += a_N;
                
            }
            b += -UdotNormal * face->get_length();
        }
        
        // * diagonal
        A_pc(c->id, c->id) += a_P;
        b_pc(c->id) += b; // * fonte
    }
    
    // + Resolvendo o problema de não ter valor prescrito de pressão...
    A_pc.row(0).zeros(); // zera a linha
    A_pc(0,0) = 1.0;
    b_pc[0] = 0.0; // prescrevendo uma correção em um nó. Deixará de ser singular a matriz...
    
    this->p_prime = arma::solve(A_pc, b_pc);
    // cout << A_pc << endl;
    // cout << b_pc << endl;
    cout << "=============================================\n";
}

/**
 * * Corrige as variáveis...
 */
void NSSolver::correct_variables(double alpha_uv, double alpha_p){
    cout << "Corrigindo as variáveis.\n";
    
    // * correção da pressão
    vector<Cell*> cells = mesh->get_cells();
    vector<Edge*> faces = mesh->get_edges();
    for(int i = 0; i < cells.size(); i++){
        // p = p^* + λ_p p'
        p[i] = p_star[i] + alpha_p * p_prime[i]; 
    }
    
    // * correção da velocidade na face
    for(int i = 0; i < faces.size(); i++){
        // * warm-up
        Edge* face = faces[i];
        pair<int,int> id_nodes = face->get_link_face_to_cell();
        int P = id_nodes.first;
        int N = id_nodes.second;
        if(!face->is_boundary_face()){
            
            double d_P = cells[P]->get_area()/a_coeff[P];
            double d_N = cells[N]->get_area()/a_coeff[N];
            double d_f = 0.5 * (d_P + d_N);

            pair<double,double> grad_p_prime_P = reconstruct_pressure_gradients(cells[P], this->p_prime);
            pair<double,double> grad_p_prime_N = reconstruct_pressure_gradients(cells[N], this->p_prime);
            
            u_face_prime[i] = d_f * 0.5 * (grad_p_prime_P.first + grad_p_prime_N.first);
            v_face_prime[i] = d_f * 0.5 * (grad_p_prime_P.second + grad_p_prime_N.second);

            // u_f = u_f^* + λ_uv u_f'
            double old_u_face = u_face[i];
            double new_u_face = old_u_face - u_face_prime[i];
            u_face[i] = old_u_face * (1 - alpha_uv) + alpha_uv * new_u_face;
            
            double old_v_face = v_face[i];
            double new_v_face = old_v_face - v_face_prime[i];
            v_face[i] = old_v_face * (1 - alpha_uv) + alpha_uv * new_v_face;
        }
    }
    
    // * correção das velocidades nos centroides
    for(int i=0; i < cells.size(); ++i){
        Cell *c = cells[i];
        double Dc = c->get_area() / a_coeff[c->id];
        pair<double,double> grad_p_prime = reconstruct_pressure_gradients(c, this->p_prime);

        // u_C = u_C^* + λ_uv u_C'
        double u_old = u_star[i];
        double u_new = u_star[i] - Dc * grad_p_prime.first;
        u[i] = u_old * (1 - alpha_uv) + alpha_uv * u_new;
        
        double v_old = v_star[i];
        double v_new = v_star[i] - Dc * grad_p_prime.second;
        v[i] = v_old * (1 - alpha_uv) + alpha_uv * v_new;
    }
    cout << "===================================================\n";


    // + atualizando
    u_star = u;
    v_star = v;
    p_star = p;

    double mass_resid = 0.0;
    for(int i = 0; i < cells.size(); i++){
        vector<int> nsigns = cells[i]->get_nsigns();
        vector<Edge*> faces = cells[i]->get_edges();
        
        double sumflux = 0;
        for(int j = 0; j < faces.size(); j++){
            Edge* f = faces[j];
            int nsign = nsigns[j];
            pair<double,double> normal = f->get_normal();
            pair<double, double> normal_corrected = make_pair(normal.first * nsign, normal.second * nsign);
        
            sumflux += (normal_corrected.first*u_face[f->id]+normal_corrected.second*v_face[f->id]) * f->get_length();
        }
        mass_resid += fabs(sumflux);
    }
    cout << "mass_resid = " << mass_resid << endl;
}

void NSSolver::export_solution(string filename, char variable){
    arma::vec var;
    if(variable == 'u')
        var = this->u_star;
    else if(variable == 'v')
        var = this->v_star;
    else if(variable == 'p')
        var = this->p_star;
    else if(variable == 'x')
        var = this->u_EXACT;
    else if(variable == 'y')
        var = this->v_EXACT;
    else if(variable == 'z')
        var = this->p_EXACT;
    
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