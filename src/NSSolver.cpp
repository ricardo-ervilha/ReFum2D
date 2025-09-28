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

    this->a = arma::vec(ncells, arma::fill::zeros);

    vector<Edge*> faces = mesh->get_edges();
    for(int i = 0; i < nfaces; i++){
        if(faces[i]->from->y == 1 && faces[i]->to->y == 1)
            u_face[i] = 1.0; // Condição do problema do lid driven cavity flow
    }
}   

NSSolver::~NSSolver(){
    // nada.
}

void NSSolver::calculate_Df(){
    cout << "Calculando os valores de Df.\n";
    vector<Edge*> faces = this->mesh->get_edges();
    for(int i = 0; i < faces.size(); ++i){
        Df[i] = mu * faces[i]->get_length() / faces[i]->get_df();
    }
    cout << "=============================================\n";
}

void NSSolver::calculate_Gf(){
    cout << "Calculando os valores de Gf.\n";
    vector<Edge*> faces = this->mesh->get_edges();
    for(int i = 0; i < faces.size(); ++i){
        Gf[i] = rho * faces[i]->get_length();
    }
    cout << "=============================================\n";
}

void NSSolver::calculate_A_mom(){
    vector<Cell*> cells = mesh->get_cells();
    for(int i = 0; i < cells.size(); ++i){
        Cell* c = cells[i];
        vector<int> &nsigns = c->get_nsigns();
        vector<Edge*> faces = c->get_edges();
        double a_P = 0;
        // contabilizando para os vizinhos
        for(int j = 0; j < faces.size(); ++j){
            int nsign = nsigns[j];
            Edge* face = faces[j];
            pair<double, double> &normal = face->get_normal();
            pair<double, double> normal_corrected = make_pair(normal.first * nsign, normal.second * nsign);
            
            int id_face = face->id;
            double U_dot_normal = normal_corrected.first * u_face[id_face] + normal_corrected.second * v_face[id_face];
            double a_N = min(0.0, -Gf[id_face] * U_dot_normal) - Df[id_face];
            if(!face->is_boundary_face()){

                a_P += max(0.0, Gf[id_face] * U_dot_normal) + Df[id_face];
                
                int id_neighbor = get_neighbor(face->get_link_face_to_cell(), c->id);
                this->A_mom(c->id, id_neighbor) += a_N;
            }else{
                a_P += max(0.0, Gf[id_face] * U_dot_normal) + Df[id_face];
            }
        }
        // contabilizando para P
        this->A_mom(c->id, c->id) += a_P;
        a[c->id] = a_P;
    }
}

void NSSolver::calculate_b_mom_x(){
    vector<Cell*> cells = mesh->get_cells();
    for(int i = 0; i < cells.size(); ++i){
        Cell* c = cells[i];
        vector<Edge*> faces = c->get_edges();
        vector<int> &nsigns = c->get_nsigns();
        double pressure_gradient = 0;
        for(int j = 0; j < faces.size(); ++j){
            int nsign = nsigns[j];
            Edge* face = faces[j];
            pair<double, double> &normal = face->get_normal();
            pair<double, double> normal_corrected = make_pair(normal.first * nsign, normal.second * nsign);
            int id_face = face->id;
            int id_neighbor = get_neighbor(face->get_link_face_to_cell(), c->id);
            
            // TODO: Interpolando com uma média aritmética.
            pressure_gradient -= ((this->p[c->id] + this->p[id_neighbor]) * 0.5) * face->get_length() * face->get_normal().first;

            if(face->is_boundary_face()){
                double U_dot_normal = normal_corrected.first * u_face[id_face] + normal_corrected.second * v_face[id_face];

                double a_N = -min(0.0, -Gf[id_face] * U_dot_normal) + Df[id_face];

                this->b_mom[c->id] += a_N * u_face[face->id];
            }
        }
        this->b_mom[c->id] += pressure_gradient + source_x * c->get_area();
    }
}

void NSSolver::calculate_b_mom_y(){
    vector<Cell*> cells = mesh->get_cells();
    for(int i = 0; i < cells.size(); ++i){
        Cell* c = cells[i];
        vector<Edge*> faces = c->get_edges();
        vector<int> &nsigns = c->get_nsigns();
        double pressure_gradient = 0;
        for(int j = 0; j < faces.size(); ++j){
            int nsign = nsigns[j];
            Edge* face = faces[j];
            pair<double, double> &normal = face->get_normal();
            pair<double, double> normal_corrected = make_pair(normal.first * nsign, normal.second * nsign);
            int id_face = face->id;
            int id_neighbor = get_neighbor(face->get_link_face_to_cell(), c->id);
            
            // TODO: Interpolando com uma média aritmética.
            pressure_gradient -= ((this->p[c->id] + this->p[id_neighbor]) * 0.5) * face->get_length() * face->get_normal().second; // * n_y

            if(face->is_boundary_face()){
                double U_dot_normal = normal_corrected.first * u_face[id_face] + normal_corrected.second * v_face[id_face];

                double a_N = -min(0.0, -Gf[id_face] * U_dot_normal) + Df[id_face];

                this->b_mom[c->id] += a_N * v_face[face->id]; // * v_face (contorno)
            }
        }
        this->b_mom[c->id] += pressure_gradient + source_y * c->get_area(); // * S_y
    }
}

void NSSolver::calculate_momentum(){

    // preenche A: serve pra u e pra v
    cout << "Calculando as equações de momento.\n";
    cout << "Calculando os coeficientes da matriz A.\n";
    calculate_A_mom();

    // preenche b para o momento em x
    cout << "Calculando os coeficientes do vetor b de x-mom.\n";
    calculate_b_mom_x();

    cout << "Calculando u_star.\n";
    this->u_star = arma::solve(this->A_mom, this->b_mom);

    // zera pra calcular o do momento em y
    this->b_mom.zeros();

    cout << "Calculando os coeficientes do vetor b de y-mom.\n";
    calculate_b_mom_y();
    cout << "Calculando v_star.\n";
    this->v_star = arma::solve(this->A_mom, this->b_mom);

    cout << "=============================================\n";
}

void NSSolver::interpolate_momentum(){
    cout << "Interpolando a velocidade nas faces.\n";
    vector<Edge*> faces = mesh->get_edges();
    vector<Cell*> cells = mesh->get_cells();
    for(int i = 0; i < faces.size(); ++i){
        Edge* face = faces[i];
        if(!face->is_boundary_face()){ // ! Não pode ser boundary face.
            pair<int,int> adjacent_cells = face->get_link_face_to_cell();
            int cell_C = adjacent_cells.first;
            int cell_F = adjacent_cells.second;

            // --- Interpolação da velocidade (U_bar) ---
            double u_avg = 0.5 * (u_star[cell_C] + u_star[cell_F]);
            double v_avg = 0.5 * (v_star[cell_C] + v_star[cell_F]);

            // --- Vetor e_CF (normalizado) ---
            pair<double,double> dCF_vec = {
                cells[cell_F]->get_centroid().first - cells[cell_C]->get_centroid().first,
                cells[cell_F]->get_centroid().second - cells[cell_C]->get_centroid().second
            };
            double dCF = face->get_df();
            pair<double,double> eCF = { dCF_vec.first / dCF, dCF_vec.second / dCF };

            // --- Distância d_f ---
            double dC = cells[cell_C]->get_area() / a[cell_C];
            double dF = cells[cell_F]->get_area() / a[cell_F];
            double d_f = 0.5 * (dC + dF);

            // --- Termo [ (pF - pC)/dCF - (gradPf ⋅ eCF) ] ---
            double pressure_diff = (p[cell_F] - p[cell_C]) / dCF;

            // gradiente de pressão interpolado na face
            auto gradP_C = reconstruct_pressure_gradients(cells[cell_C], this->p);
            auto gradP_F = reconstruct_pressure_gradients(cells[cell_F], this->p);
            pair<double,double> gradP_face = {
                0.5 * (gradP_C.first + gradP_F.first),
                0.5 * (gradP_C.second + gradP_F.second)
            };

            // produto escalar gradPf ⋅ eCF
            double grad_dot_eCF = gradP_face.first * eCF.first + gradP_face.second * eCF.second;

            double correction = (pressure_diff - grad_dot_eCF);

            // --- Correção Rhie–Chow ---
            double u_corr = d_f * correction * eCF.first;
            double v_corr = d_f * correction * eCF.second;

            u_face[face->id] = u_avg - u_corr;
            v_face[face->id] = v_avg - v_corr;
        }
    }
    cout << "=============================================\n";
}

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
            // TODO: Esse tratamento faz sentido ?
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

void NSSolver::pressure_correction_poisson(){
    cout << "Resolvendo a poisson para encontrar a correção da pressão.\n";
    vector<Cell*> cells = mesh->get_cells();
    for(int i = 0; i < cells.size(); ++i){
        Cell* c = cells[i];
        vector<Edge*> faces = c->get_edges();
        vector<int> nsigns = c->get_nsigns();
        double a_P = 0;
        double b = 0;
        
        for(int j = 0; j < faces.size(); ++j){
            Edge* face = faces[j];
            int nsign = nsigns[j];
            pair<double,double> normal = face->get_normal();
            pair<double, double> normal_corrected = make_pair(normal.first * nsign, normal.second * nsign);
            
            if(!face->is_boundary_face()){
                int id_neighbor = get_neighbor(face->get_link_face_to_cell(), c->id);
            
                double U_dot_normal = normal_corrected.first * u_face[face->id]  + normal_corrected.second * v_face[face->id];
                // double a_N = min(0.0, -Gf[face->id] * U_dot_normal) - Df[face->id];
                
                double d_P = cells[c->id]->get_area()/a[c->id];
                double d_N = cells[id_neighbor]->get_area()/a[id_neighbor];
                double d_f = 0.5 * (d_P + d_N);

                double dcfx = cells[id_neighbor]->get_centroid().first - cells[c->id]->get_centroid().first;
                double dcfy = cells[id_neighbor]->get_centroid().second - cells[c->id]->get_centroid().second;

                double sfx = face->get_length() * normal_corrected.first;
                double sfy = face->get_length() * normal_corrected.second;

                double coeff = (dcfx * d_f * sfx + dcfy * d_f * sfy)/(dcfx*dcfx + dcfy*dcfy);
                
                this->A_pc(c->id, id_neighbor) += -coeff;
                a_P += coeff;
                b += U_dot_normal * face->get_length();
            }else{
                double U_dot_normal = normal_corrected.first * u_face[face->id]  + normal_corrected.second * v_face[face->id];
                this->b_pc(c->id) += - U_dot_normal * face->get_length();
            }
        }
        /* Resolvendo o problema de não ter valor prescrito de pressão...*/
        A_pc.row(0).zeros(); // zera a linha
        A_pc(0,0) = 1.0;
        b_pc[0] = 0.0; // prescrevendo uma correção em um nó. Deixará de ser singular a matriz...

        this->A_pc(c->id, c->id) = a_P;
        this->b_pc(c->id) += -b;
    }
    this->p_prime = arma::solve(A_pc, b_pc);
    cout << "=============================================\n";
}

void NSSolver::pressure_correction_poisson_modified(){
    cout << "Resolvendo a poisson para encontrar a correção da pressão.\n";
    vector<Cell*> cells = mesh->get_cells();
    for(int i = 0; i < cells.size(); ++i){
        Cell* c = cells[i];
        vector<Edge*> faces = c->get_edges();
        vector<int> nsigns = c->get_nsigns();
        
        double a_P = 0;
        double b = 0;
        
        for(int j = 0; j < faces.size(); ++j){
            Edge* face = faces[j];
            int nsign = nsigns[j];
            pair<double,double> normal = face->get_normal();
            pair<double, double> normal_corrected = make_pair(normal.first * nsign, normal.second * nsign);

            if(!face->is_boundary_face()){
                int id_neighbor = get_neighbor(face->get_link_face_to_cell(), c->id);
                double U_dot_normal = normal_corrected.first * u_face[face->id]  + normal_corrected.second * v_face[face->id];

                double d_N = cells[id_neighbor]->get_area() / a[id_neighbor];
                double d_P = cells[c->id]->get_area() / a[c->id];
                double d_f = 0.5* (d_N + d_P);

                double a_N = d_f * face->get_length();
                this->A_pc(c->id, id_neighbor) += -a_N;
                a_P += a_N;

                b += U_dot_normal * face->get_length();
            }else{
                // P'_b = P'_c.
                int id_neighbor = get_neighbor(face->get_link_face_to_cell(), c->id);
                double U_dot_normal = normal_corrected.first * u_face[face->id]  + normal_corrected.second * v_face[face->id];

                double d_P = cells[c->id]->get_area() / a[c->id];

                double a_N = d_P * face->get_length();
                this->A_pc(c->id, c->id) += -a_N;
                a_P += a_N;

                b += U_dot_normal * face->get_length();
            }
        }
        A_pc(c->id, c->id) += a_P;
        b_pc(c->id) += b;
    }
    cout << this->A_pc << endl;
    this->p_prime = arma::solve(A_pc, b_pc);
    cout << "=============================================\n";
}

void NSSolver::correct_variables(){
    cout << "Corrigindo as variáveis.\n";
    vector<Cell*> cells = mesh->get_cells();
    vector<Edge*> faces = mesh->get_edges();
    for(int i = 0; i < cells.size(); i++){
        p[i] = p_star[i] + p_prime[i]; // correto.
    }
    
    for(int i = 0; i < faces.size(); i++){
        Edge* face = faces[i];
        pair<int,int> id_nodes = face->get_link_face_to_cell();
        int P = id_nodes.first;
        int N = id_nodes.second;
        if(!face->is_boundary_face()){
            
            double d_P = cells[P]->get_area()/a[P];
            double d_N = cells[N]->get_area()/a[N];
            double d_f = 0.5 * (d_P + d_N);
            
            u_face_prime[i] = d_f * (p_prime[P] - p_prime[id_nodes.second])/face->get_df();
            v_face_prime[i] = d_f * (p_prime[P] - p_prime[id_nodes.second])/face->get_df();
        }
    }
    
    for(int i = 0; i < faces.size(); i++){
        u_face[i] = u_face[i] - u_face_prime[i] * faces[i]->get_length() * faces[i]->get_normal().first;
        v_face[i] = v_face[i] - v_face_prime[i] * faces[i]->get_length() * faces[i]->get_normal().second;
    }
    
    for(int i=0; i < cells.size(); ++i){
        Cell *c = cells[i];
        double Dc = c->get_area() / a[c->id];
        pair<double,double> grad_p_prime = reconstruct_pressure_gradients(c, this->p_prime);

        u_star[i] = u_star[i] - Dc * grad_p_prime.first;
        v_star[i] = v_star[i] - Dc * grad_p_prime.second;
    }
    cout << "===================================================\n";
}

void NSSolver::export_solution(string filename, char variable){
    arma::vec var;
    if(variable == 'u')
        var = this->u_star;
    else if(variable == 'v')
        var = this->v_star;
    else if(variable == 'p')
        var = this->p;
    
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