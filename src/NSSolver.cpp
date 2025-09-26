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
    this->A_pc = arma::sp_mat(ncells, ncells);
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

    cout << this->A_mom << endl;

    // preenche b para o momento em x
    cout << "Calculando os coeficientes do vetor b de x-mom.\n";
    calculate_b_mom_x();

    cout << this->b_mom << endl;

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
    vector<Edge*> faces = mesh->get_edges();
    vector<Cell*> cells = mesh->get_cells();
    for(int i = 0; i < faces.size(); ++i){
        Edge* face = faces[i];
        if(!face->is_boundary_face()){ // ! Não pode ser boundary face.
            pair<int,int> id_nodes = face->get_link_face_to_cell();
            int P = id_nodes.first;
            int N = id_nodes.second;
            
            // TODO: Interpolando com média aritmética.
            double u_bar = 0.5 * (u_star[P] + u_star[N]);
            double v_bar = 0.5 * (v_star[P] + v_star[N]);
            
            double U_dot_normal = face->get_normal().first * u_face[face->id]  + face->get_normal().second * v_face[face->id];
            double a_N = min(0.0, -Gf[face->id] * U_dot_normal) - Df[face->id];

            double d_P = cells[P]->get_area()/a[P];
            double d_N = cells[N]->get_area()/a_N;
            double d_f = 0.5 * (d_P + d_N);

            double dx = fabs(face->from->x - face->to->x);
            double dy = fabs(face->from->y - face->to->y);
            double partial_p_partial_x_f = (this->p[N] - this->p[P])/dx;
            double partial_p_partial_y_f = (this->p[N] - this->p[P])/dy;

            pair<double,double> gradient_pressure_P = reconstruct_pressure_gradients(cells[P]);
            pair<double,double> gradient_pressure_N = reconstruct_pressure_gradients(cells[N]);

            u_face[i] = u_bar - d_f * (partial_p_partial_x_f + 0.5 * (gradient_pressure_P.first + gradient_pressure_N.first));
            v_face[i] = v_bar - d_f * (partial_p_partial_y_f + 0.5 * (gradient_pressure_P.second + gradient_pressure_N.second));
        }
    }
}

pair<double,double> NSSolver::reconstruct_pressure_gradients(Cell *c){  

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
            y[j] = this->p[nb] - this->p[c->id];  
        }
    }
        
    // + Resolve sistema sobre-determinado M x = y.
    arma::vec x = arma::solve(M, y);
    double px = x[0];
    double py = x[1];
    
    return make_pair(px, py);
}

void NSSolver::pressure_correction_poisson(){
    vector<Cell*> cells = mesh->get_cells();
    for(int i = 0; i < cells.size(); ++i){
        Cell* c = cells[i];
        vector<Edge*> faces = c->get_edges();
        double a_P = 0;
        double b = 0;
        for(int j = 0; j < faces.size(); ++j){
            Edge* face = faces[i];
            int id_neighbor = get_neighbor(face->get_link_face_to_cell(), c->id);
            
            double U_dot_normal = face->get_normal().first * u_face[face->id]  + face->get_normal().second * v_face[face->id];
            double a_N = min(0.0, -Gf[face->id] * U_dot_normal) - Df[face->id];
            
            double d_P = cells[c->id]->get_area()/a[c->id];
            double d_N = cells[id_neighbor]->get_area()/a_N;
            double d_f = 0.5 * (d_P + d_N);

            this->A_pc(c->id, id_neighbor) += - d_f * face->get_length();
            a_P += d_f * face->get_length();

            b += U_dot_normal * face->get_length();
        }
        this->A_pc(c->id, c->id) += a_P;
        this->b_pc(c->id) += - b;
    }
    this->p_prime = arma::spsolve(A_pc, b_pc);
}

void NSSolver::correct_variables(){
    vector<Cell*> cells = mesh->get_cells();
    vector<Edge*> faces = mesh->get_edges();
    for(int i = 0; i < cells.size(); i++){
        p[i] = p_star[i] + p_prime[i];
    }

    for(int i = 0; i < faces.size(); i++){
        Edge* face = faces[i];
        pair<int,int> id_nodes = face->get_link_face_to_cell();
        int P = id_nodes.first;

        int id_neighbor = get_neighbor(face->get_link_face_to_cell(), P);
        double U_dot_normal = face->get_normal().first * u_face[face->id]  + face->get_normal().second * v_face[face->id];
        double a_N = min(0.0, -Gf[face->id] * U_dot_normal) - Df[face->id];
        
        double d_P = cells[P]->get_area()/a[P];
        double d_N = cells[id_neighbor]->get_area()/a_N;
        double d_f = 0.5 * (d_P + d_N);
        
        u_face_prime[i] = d_f * (p_prime[P] - p_prime[id_nodes.second]);
        v_face_prime[i] = d_f * (p_prime[P] - p_prime[id_nodes.second]);
    }

    for(int i = 0; i < faces.size(); i++){
        u_face[i] = u_face[i] + u_face_prime[i];
        v_face[i] = v_face[i] + v_face_prime[i];
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
    vtk_file << "SCALARS Temperatura double 1\n";
    vtk_file << "LOOKUP_TABLE default\n";
    for(int i = 0; i < cells.size(); i++){
        vtk_file << this->v_star[i] << endl;
    }

    vtk_file.close();
}