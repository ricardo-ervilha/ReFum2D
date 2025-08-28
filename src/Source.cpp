#include "Source.h"
#include "FVMSolver.h"

Source::Source(FVMSolver* s, double (*source)(double, double)){
    this->solver = s;
    this->sourcefunc = source;
    this->sourcecc = arma::vec(solver->mesh->get_ncells(), arma::fill::zeros);
    
    this->interpolate_sources();
}

Source::~Source(){
    //Nada
}

void Source::interpolate_sources(){
    vector<Cell*> cells = solver->mesh->get_cells();

    for(Cell* c : cells){
        vector<Node*> nodes = c->get_nodes();
        double num = 0;
        double den = 0;
        pair<double,double> centroidP = c->get_centroid();
        for(Node* n: nodes){
            // + Calcula termo fonte para o nó da malha que é vertice do VC
            double source_node = this->sourcefunc(n->x, n->y);
            
            // + Calcula distancia entre nó e centro da celula
            double dist_center_cell_to_node = distance(
                centroidP.first, centroidP.second,
                n->x, n->y
            );

            // + Em cima vai ponderado o termo pela distancia. Em baixo está a soma das distancias
            num += source_node * dist_center_cell_to_node;
            den += dist_center_cell_to_node;
        }

        this->sourcecc[c->id] = num / den;
    }
}

void Source::assemblyCoefficients(){
    vector<Cell*> cells = solver->mesh->get_cells();

    for(Cell* c : cells){
        // TODO: sinal
        solver->b[c->id] += - this->sourcecc[c->id] * c->get_area();
    }
}