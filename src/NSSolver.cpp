#include "NSSolver.h"
#include "Mesh.h"
#include "Edge.h"
#include "Cell.h"

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