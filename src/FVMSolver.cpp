#include "../include/FVMSolver.h"
#include "GradientReconstruction.h"

FVMSolver::FVMSolver(Mesh *mesh, BoundaryCondition *bc1, BoundaryCondition* bc2, BoundaryCondition* bc3, BoundaryCondition *bc4){
    // * Salva a malha que vêm pré-processada da main.
    this->mesh = mesh;

    // * Salva as condições de contorno em ordem anti-horária: bottom, right, top, left.
    this->boundaries.emplace(bc1->get_location(), bc1);
    this->boundaries.emplace(bc2->get_location(), bc2);
    this->boundaries.emplace(bc3->get_location(), bc3);
    this->boundaries.emplace(bc4->get_location(), bc4);

    BoundaryCondition::apply_mins_maxs(mesh->get_xmin(), mesh->get_xmax(), mesh->get_ymin(), mesh->get_ymax());

    // ========================================================================================
    // ========================================================================================
    // ! Inicialização das estruturas de dados do problema.
    int ncells = this->mesh->get_ncells();

    this->A = arma::sp_mat(ncells, ncells);
    this->b = arma::vec(ncells, arma::fill::zeros);
    this->b_aux = arma::vec(ncells, arma::fill::zeros);
    this->u_old = arma::vec(ncells, arma::fill::zeros);  
    this->u_new = arma::vec(ncells, arma::fill::zeros);  
    this->gradients = arma::mat(ncells, 2);
}

FVMSolver::~FVMSolver(){
    //Nada
}

void FVMSolver::SteadySolver(Diffusion* d, bool cd, int num_iter_cd){
    GradientReconstruction* g = new GradientReconstruction(this);
    if(!cd)
    {
        this->u_new = arma::spsolve(A, b);
    }    
    else{
        for(int i = 0; i < num_iter_cd; i++){
            b_aux.zeros(); // ZERA o explícito
            g->reconstruct_gradients(); // reconstroi os gradientes
            d->cross_diffusion(); // Computa a difusão cruzada;
            b_aux = b_aux + b; // soma o vetor com contribuição fixa com o que possui a contribuição corrigida
            this->u_new = arma::spsolve(A, b_aux); // resolve com a correção  
        }
    }
}

/**
 * * Computa erro usando norma L_2.
 */
void FVMSolver::compute_error(double (*exact)(double, double)){
    arma::vec exact_vect(mesh->get_ncells());
    vector<Cell*> cells = mesh->get_cells();
    for(int i = 0; i < cells.size(); i++){ //para cada célula
        pair<double,double>& centroid = cells[i]->get_centroid(); // obtém o centroide da celula
        exact_vect[i] = exact(centroid.first, centroid.second); // ! Aplica a exata no centroide 
    }

    double norml2 = 0;
    double sumAreas = 0;
    for (int i = 0; i < cells.size(); ++i) {
        norml2 += (u_new[i] - exact_vect[i]) * (u_new[i] - exact_vect[i]) * cells[i]->get_area();
        sumAreas +=  cells[i]->get_area();
    }
    norml2 = sqrt(norml2/sumAreas);
    
    cout << "\nNorma L2: " << norml2 << endl;
}

void FVMSolver::export_solution(string filename){
    
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
        vtk_file << this->u_new[i] << endl;
    }

    vtk_file.close();
}