#include "Diffusion.h"
#include "FVMSolver.h"
#include "Mesh.h"
#include "Cell.h"

Diffusion::Diffusion(FVMSolver* s, double (*g)(double, double)){
    this->solver = s;
    this->gammafunc = g;
    this->gammaf = arma::vec(solver->mesh->get_nedges(), arma::fill::zeros);
    
    this->interpolate_gammas(); // calcula gamma nas faces
}
 

Diffusion::~Diffusion(){
    //Nada
}

void Diffusion::interpolate_gammas(){
    vector<Edge*> edges = solver->mesh->get_edges();
    for(Edge *e : edges){
        const Node* from = e->from; 
        const Node* to = e->to;
        
        double gamma_from = gammafunc(from->x, from->y);
        double gamma_to = gammafunc(to->x, to->y);

        // valor de gamma no centro da face interpolado por uma média aritmética
        this->gammaf[e->id] = 0.5 * (gamma_from + gamma_to);
    }
}

/**
 * * Função para descobrir o vetor n1
 * @param: p é a celula que está sendo avaliada
 * @param: n é a celula vizinha
 * @param: normal é a normal da face apontando para fora da célula
 */
pair<double,double> Diffusion::compute_n1(pair<double,double>& p, pair<double,double>& n, pair<double,double>& normal){
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


void Diffusion::assembleCoefficients(){
    vector<Cell*>& cells = solver->mesh->get_cells();
    for(auto cell: cells)
        this->control_volume_non_orthogonal_diffusion(cells, cell);
}

void Diffusion::control_volume_non_orthogonal_diffusion(vector<Cell*> cells, Cell *c){
    // * Faces da célula e sinais da normal para apontar para fora corretamente.
    vector<Edge*> faces = c->get_edges();
    vector<int> &nsigns = c->get_nsigns();

    // * Centroide da célula que está sendo avaliada
    pair<double,double>& centroidP = c->get_centroid();
    
    for (int i = 0; i < faces.size(); ++i) {
        Edge* face = faces[i];
        int nsign = nsigns[i];

        pair<double, double> &middleFace = face->get_middle();
        pair<double, double> &normal = face->get_normal();
        pair<double, double> normal_corrected = make_pair(normal.first * nsign, normal.second * nsign);

        if(face->is_boundary_face()){
            BoundaryLocation local = BoundaryCondition::find_location(face);
            BoundaryType bt = solver->boundaries[local]->get_type();
            double bound_value = solver->boundaries[local]->apply(middleFace.first, middleFace.second);
            
            if(bt == DIRICHLET){
                // TODO: Rever com o professor essa conta...

                // * Quando é dirichlet N1 ligará centro de P e centro da face.
                pair<double,double> n1 = compute_n1(centroidP, middleFace, normal_corrected);
                double normn1 = sqrt(n1.first * n1.first + n1.second * n1.second); // |n1|

                solver->A(c->id, c->id) += (gammaf[face->id] * face->get_length() * normn1) / face->get_df();

                solver->b[c->id] += (gammaf[face->id] * face->get_length() * normn1 * bound_value) / face->get_df();

            }else if(bt == NEUMANN){
                // * Nesse caso é mais simples pois já tem o valor daquele gradiente produto com normal.
                solver->b(c->id) += (gammaf[face->id] * face->get_length() * bound_value);
            }
        }else{
            int nb = get_neighbor(face->get_link_face_to_cell(), c->id);
            pair<double,double>& centroidN = cells[nb]->get_centroid(); // centroide do vizinho
            
            // + Vetor n1 unitário que liga o centroide de P com N.
            pair<double,double> n1 = compute_n1(centroidP, centroidN, normal_corrected);

            double norm_n1 = sqrt(n1.first * n1.first + n1.second * n1.second); // + |n1|

            // contribuição na diagonal (phi_P)
            solver->A(c->id, c->id) += (gammaf[face->id] * face->get_length() * norm_n1) / face->get_df();

            // contribuição no vizinho (phi_N)
            solver->A(c->id, nb) += - (gammaf[face->id] * face->get_length() * norm_n1) / face->get_df();
        }
    }
}

void Diffusion::cross_diffusion(){
    vector<Cell*> cells = solver->mesh->get_cells();

    for(int i = 0; i < cells.size(); i++){ 
        vector<Edge*> facesOfCell = cells[i]->get_edges();
        vector<int>& nsigns = cells[i]->get_nsigns();
        pair<double,double>& centroidP = cells[i]->get_centroid();

        for(int j = 0; j < facesOfCell.size(); j++){
            Edge* e = facesOfCell[j];
            
            int nsign = nsigns[j];

            pair<double, double> &middleFace = e->get_middle();
            pair<double, double>& normal = e->get_normal();
            pair<double,double> normal_corrected = make_pair(normal.first * nsign, normal.second * nsign);    

            if(e->is_boundary_face()){
                BoundaryLocation local = BoundaryCondition::find_location(e);
                BoundaryType bt = solver->boundaries[local]->get_type();
                double bound_value = solver->boundaries[local]->apply(middleFace.first, middleFace.second);
                
                if(bt == DIRICHLET){
                    pair<double,double> n1 = compute_n1(centroidP, middleFace, normal_corrected);

                    // + encontrando o n2: n1 + n2 = nf => n2 = nf - n1
                    double n2x = normal_corrected.first - n1.first;
                    double n2y = normal_corrected.second - n1.second;

                    // estou assumindo que o gradiente da face é o gradiente da celula.
                    double graddotnormal = solver->gradients(i,0) * n2x + solver->gradients(i,1) * n2y;

                    solver->b_aux[i] += (gammaf[e->id] * e->get_length() * graddotnormal);
                }else if(bt == NEUMANN){
                    // * não precisa tratar
                }
            }else{
                int nb = get_neighbor(e->get_link_face_to_cell(), cells[i]->id);
                pair<double,double>& centroidN = cells[nb]->get_centroid(); // centroide do vizinho

                pair<double,double> n1 = compute_n1(centroidP, centroidN, normal_corrected);

                double n2x = normal_corrected.first - n1.first;
                double n2y = normal_corrected.second - n1.second;

                /*Calculo da distancia de P até center(gface) & N até center(gface)*/
                pair<double,double>& midFace = facesOfCell[j]->get_middle();
                double d1 = distance(centroidP.first, centroidP.second, midFace.first, midFace.second);
                double d2 = distance(centroidN.first, centroidN.second, midFace.first, midFace.second);
                
                // + Interpolação do gradiente da face
                double gradx = (solver->gradients(i,0)*d2 + solver->gradients(nb,0)*d1)/(d1+d2);
                double grady = (solver->gradients(i,1)*d2 + solver->gradients(nb,1)*d1)/(d1+d2);

                double graddotnormal = gradx * n2x + grady * n2y;

                solver->b_aux[i] += (gammaf[e->id] * facesOfCell[j]->get_length() * graddotnormal);
            }
        }
    }
}