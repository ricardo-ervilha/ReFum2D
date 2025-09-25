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
    this->A_mom = arma::sp_mat(ncells, ncells);
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

    // variáveis relacionadas a pressão
    this->p = arma::vec(ncells, arma::fill::zeros);
    this->p_prime = arma::vec(ncells, arma::fill::zeros);
    this->p_star = arma::vec(ncells, arma::fill::zeros);

    // coeficientes de Difusão e Convecção
    this->Df = arma::vec(nfaces, arma::fill::zeros);
    this->Gf = arma::vec(nfaces, arma::fill::zeros);

    this->a = arma::vec(ncells, arma::fill::zeros);
}   

NSSolver::~NSSolver(){
    // nada.
}

void NSSolver::calculate_Df(){
    vector<Edge*> faces = this->mesh->get_edges();
    for(int i = 0; i < faces.size(); ++i){
        Df[i] = mu * faces[i]->get_length() / faces[i]->get_df();
    }
}

void NSSolver::calculate_Gf(){
    vector<Edge*> faces = this->mesh->get_edges();
    for(int i = 0; i < faces.size(); ++i){
        // TODO: essa normal da face pode ter dois sentidos. Qual deveria considerar ?
        double Uf_dot_nf = u_face[i] * faces[i]->get_normal().first + u_face[i] * faces[i]->get_normal().second;
        Gf[i] = rho * faces[i]->get_length() * Uf_dot_nf;
    }
}

void NSSolver::calculate_A_mom(){
    vector<Cell*> cells = mesh->get_cells();
    for(int i = 0; i < cells.size(); ++i){
        Cell* c = cells[i];
        vector<Edge*> faces = c->get_edges();
        double a_P = 0;
        // contabilizando para os vizinhos
        for(int j = 0; j < faces.size(); ++j){
            Edge* face = faces[j];
            int id_face = face->id;

            double a_N = min(0.0, -Gf[id_face]) - Df[id_face];
            a_P += max(0.0, Gf[id_face]) + Df[id_face];
            
            int id_neighbor = get_neighbor(face->get_link_face_to_cell(), c->id);
            this->A_mom(c->id, id_neighbor) += a_N;
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
        
        double pressure_gradient = 0;
        for(int j = 0; j < faces.size(); ++j){
            Edge* face = faces[j];
            int id_face = face->id;
            int id_neighbor = get_neighbor(face->get_link_face_to_cell(), c->id);
            // TODO: Interpolando com uma média aritmética.
            pressure_gradient -= ((this->p[c->id] + this->p[id_neighbor]) * 0.5) * face->get_length() * face->get_normal().first;
        }
        this->b_mom[c->id] = pressure_gradient + source_x * c->get_area();
    }
}

void NSSolver::calculate_b_mom_y(){
    vector<Cell*> cells = mesh->get_cells();
    for(int i = 0; i < cells.size(); ++i){
        Cell* c = cells[i];
        vector<Edge*> faces = c->get_edges();
        
        double pressure_gradient = 0;
        for(int j = 0; j < faces.size(); ++j){
            Edge* face = faces[j];
            int id_face = face->id;
            int id_neighbor = get_neighbor(face->get_link_face_to_cell(), c->id);
            // TODO: Interpolando com uma média aritmética.
            pressure_gradient -= ((this->p[c->id] + this->p[id_neighbor]) * 0.5) * face->get_length() * face->get_normal().second;
        }
        this->b_mom[c->id] = pressure_gradient + source_y * c->get_area();
    }
}

void NSSolver::calculate_momentum(){

    // preenche A: serve pra u e pra v
    calculate_A_mom();

    // preenche b para o momento em x
    calculate_b_mom_x();

    this->u_star = arma::spsolve(this->A_mom, this->b_mom);

    // zera pra calcular o do momento em y
    this->b_mom.zeros();

    calculate_b_mom_y();
    this->v_star = arma::spsolve(this->A_mom, this->b_mom);
}

void NSSolver::interpolate_momentum(){
    vector<Edge*> faces = mesh->get_edges();
    vector<Cell*> cells = mesh->get_cells();
    for(int i = 0; i < faces.size(); ++i){
        Edge* face = faces[i];
        pair<int,int> id_nodes = face->get_link_face_to_cell();
        int P = id_nodes.first;
        int N = id_nodes.second;
        
        // TODO: Interpolando com média aritmética.
        double u_bar = 0.5 * (u_star[P] + u_star[N]);
        double v_bar = 0.5 * (v_star[P] + v_star[N]);
        
        double d_P = cells[P]->get_area()/a[P];
        double d_N = cells[N]->get_area()/a[N];
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

            double d_P = cells[c->id]->get_area()/a[c->id];
            double d_N = cells[id_neighbor]->get_area()/a[id_neighbor];
            double d_f = 0.5 * (d_P + d_N);

            this->A_pc(c->id, id_neighbor) += - d_f * face->get_length();
            a_P += d_f * face->get_length();

            double U_dot_normal = u_face[face->id] * face->get_normal().first + v_face[face->id] * face->get_normal().second;
            b += U_dot_normal * face->get_length();
        }
        this->A_pc(c->id, c->id) += a_P;
        this->b_pc(c->id) += - b;
    }
    this->p_prime = arma::spsolve(A_pc, b_pc);
}