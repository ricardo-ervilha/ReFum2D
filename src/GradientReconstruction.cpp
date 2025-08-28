#include "GradientReconstruction.h"
#include "FVMSolver.h"

GradientReconstruction::GradientReconstruction(FVMSolver *s){
    this->solver = s;
}

GradientReconstruction::~GradientReconstruction(){
    //Nada
}

void GradientReconstruction::reconstruct_gradients(){
    vector<Cell*>& cells = solver->mesh->get_cells();
    
    for(int i = 0; i < cells.size(); i++){ 

        Cell *c = cells[i];

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
                BoundaryLocation local = BoundaryCondition::find_location(e);
                BoundaryType bt = solver->boundaries[local]->get_type();
                double bound_value = solver->boundaries[local]->apply(middleFace.first, middleFace.second);
            
                if(bt == DIRICHLET){
                    double dx = middleFace.first - centroidP.first;
                    double dy = middleFace.second - centroidP.second;
                    // + [dx dy]
                    M(j,0) = dx;
                    M(j, 1) = dy;

                     // + [u_B - u_P]
                    BoundaryLocation local = BoundaryCondition::find_location(edgesOfCell[j]);
                    y[j] = bound_value - solver->u_new[i];  

                }else if(bt == NEUMANN){
                    // * No caso de neumann já tem o gradiente...
                }
            }else{
                int nb = get_neighbor(e->get_link_face_to_cell(), c->id);
                pair<double,double>& centroidN = cells[nb]->get_centroid();

                double dx = centroidN.first - centroidP.first;
                double dy = centroidN.second - centroidP.second;
                
                // + [dx dy]
                M(j,0) = dx;
                M(j, 1) = dy;

                // + [u_N - u_P]
                y[j] = solver->u_new[nb] - solver->u_new[i];  
            }
        }
        
        // + Resolve sistema sobre-determinado M x = y.
        arma::vec x = arma::solve(M, y);
        this->solver->gradients(i,0) = x[0];
        this->solver->gradients(i,1) = x[1];
    }
}