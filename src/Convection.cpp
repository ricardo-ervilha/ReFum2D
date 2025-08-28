#include "Convection.h"
#include "FVMSolver.h"

Convection::Convection(FVMSolver *s, double (*r)(double, double), pair<double, double> (*ufunc)(double, double)){
    this->solver = s;
    this->rhofunc = r;
    this->Ufunc = ufunc;

    this->rhof = arma::vec(solver->mesh->get_nedges(), arma::fill::zeros);
    this->Uf = arma::mat(solver->mesh->get_nedges(), 2);
    
    this->interpolate_properties(); // calcula gamma nas faces
}

Convection::~Convection(){
    //Nada
}

void Convection::interpolate_properties(){
    vector<Edge*> edges = solver->mesh->get_edges();
    for(Edge *e : edges){
        const Node* from = e->from; 
        const Node* to = e->to;
        
        double rho_from = rhofunc(from->x, from->y);
        double rho_to = rhofunc(to->x, to->y);

        pair<double,double> U_from = Ufunc(from->x, from->y);
        pair<double,double> U_to = Ufunc(to->x, to->y);

        // valor de gamma no centro da face interpolado por uma média aritmética
        this->rhof[e->id] = 0.5 * (rho_from + rho_to);

        // valor da velocidade x interpolado
        this->Uf(e->id, 0) = 0.5 * (U_from.first + U_to.first);

        // valor da velocidade y interpolado
        this->Uf(e->id, 1) = 0.5 * (U_from.second + U_to.second);
    }
}

void Convection::assembleCoefficients(){
    vector<Cell*>& cells = solver->mesh->get_cells();
    for(int i = 0; i < cells.size(); i++){
        this->control_volume_convection(cells, cells[i]);
    }
}

void Convection::control_volume_convection(vector<Cell*> cells, Cell* P){
    
    vector<Edge *> faces = P->get_edges(); // * faces/lados da célula
    vector<int> &nsigns = P->get_nsigns();
    
    for (int i = 0; i < faces.size(); i++)
    {
        Edge *face = faces[i];
        int nsign = nsigns[i]; // sinal corrigida da normal

        // centro da face & normal da face corrigida
        pair<double, double> &middleFace = face->get_middle();
        pair<double, double> &normal = face->get_normal();
        pair<double, double> normal_corrected = make_pair(normal.first * nsign, normal.second * nsign);

        double UdotNormal = Uf(face->id, 0) * normal_corrected.first + Uf(face->id, 1) * normal_corrected.second;
        double G_f = rhof[face->id] * UdotNormal * face->get_length();
        
        if (face->is_boundary_face())
        {
            if (G_f > 0)
            {
                solver->A(P->id, P->id) += G_f;
            }
            else
            {
                BoundaryLocation local = BoundaryCondition::find_location(face);
                BoundaryType bt = solver->boundaries[local]->get_type();
                double bound_value = solver->boundaries[local]->apply(middleFace.first, middleFace.second);
            
                if (bt == DIRICHLET)
                {
                    solver->b[P->id] += -G_f * bound_value;
                }
                else if (bt == NEUMANN)
                {
                    // + usar série de Taylor
                    // * phi_f = phi_p + (\partial phi/ \partial n) \delta n
                    // * ficará então: Gf * [phi_p + (\partial phi/ \partial n) \delta n]
                    // ! = Gf * phi_p + Gf * \partial phi/ \partial n) \delta n
                    // ? Gf * phi_p entra na A
                    // ? Gf * \partial phi/ \partial n) \delta n entra na b

                    solver->A(P->id, P->id) += G_f;

                    solver->b[P->id] += -G_f * bound_value * face->get_df();
                }
            }
        }
        else
        {
            if (G_f > 0)
            {
                // +phiP
                solver->A(P->id, P->id) += G_f;
            }
            else
            {
                // +phiNeighbor
                int nb = get_neighbor(face->get_link_face_to_cell(), P->id);
                solver->A(P->id, nb) += G_f;
            }
        }
    }
}

void Convection::linear_upwind(){
    // TODO: implementar o linear_upwind usando o b_aux e a reconstrução do gradiente e verificar o erro.
}