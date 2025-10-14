#include "NSSolver.h"
#include "Mesh.h"
#include "Edge.h"
#include "Cell.h"
#include "BoundaryCondition.h"
#include "utils.h"

NSSolver::NSSolver(Mesh *mesh, float mu, float rho){
    PetscInitialize(NULL, NULL, NULL, NULL);
    
    KSPCreate(PETSC_COMM_WORLD, &ksp);

    this->mesh = mesh;
    this->mu = mu;
    this->rho = rho;

    // ! Inicialização de todas as variáveis que NSSolver guarda.
    int ncells = this->mesh->get_ncells();
    int nfaces = this->mesh->get_nedges();
    
    //Criação dos vetores do PETSC
    VecCreate(PETSC_COMM_WORLD, &b_mom_x);
    VecSetSizes(b_mom_x, PETSC_DECIDE, ncells);
    VecSetFromOptions(b_mom_x); // DEFINE b_mom_x como vetor padrão contendo n_cells
    
    // criando vetores com mesmo layout de b_mom_x
    VecDuplicate(b_mom_x, &b_mom_y);
    VecDuplicate(b_mom_x, &b_pc);
    VecDuplicate(b_mom_x, &uc);
    VecDuplicate(b_mom_x, &vc);
    VecDuplicate(b_mom_x, &pc);
    VecDuplicate(b_mom_x, &uc_aux);
    VecDuplicate(b_mom_x, &vc_aux);
    VecDuplicate(b_mom_x, &pc_aux);
    VecDuplicate(b_mom_x, &ap);
    VecDuplicate(b_mom_x, &pcorr);
    VecDuplicate(b_mom_x, &ucorr);
    VecDuplicate(b_mom_x, &vcorr);

    VecCreate(PETSC_COMM_WORLD, &u_face);
    VecSetSizes(u_face, PETSC_DECIDE, nfaces);
    VecSetFromOptions(u_face); // DEFINE u_face como vetor padrão contendo n_faces
    VecDuplicate(u_face, &v_face);
    VecDuplicate(u_face, &p_face);
    VecDuplicate(u_face, &mdotf);
    VecDuplicate(u_face, &mdotfcorr);
    VecDuplicate(u_face, &pfcorr);
    VecDuplicate(u_face, &wf);

    MatCreate(PETSC_COMM_WORLD, &A_mom);
    MatSetSizes(A_mom, PETSC_DECIDE, PETSC_DECIDE, ncells, ncells);
    MatSetFromOptions(A_mom);
    MatSetUp(A_mom);

    MatCreate(PETSC_COMM_WORLD, &A_pc);
    MatSetSizes(A_pc, PETSC_DECIDE, PETSC_DECIDE, ncells, ncells);
    MatSetFromOptions(A_pc);
    MatSetUp(A_pc);

    // =======================================================================================
    // TODO: Preenche boundary top do lid-driven cavity flow.
    vector<Edge*> faces = mesh->get_edges();
    for(PetscInt i = 0; i < nfaces; i++){
        Edge* face = faces[i];
        if(face->from->y == 1.0 && face->to->y == 1.0){
            // Top 
            VecSetValue(u_face, face->id, 1.0, INSERT_VALUES);
        }
    }

    this->compute_wf();
}   

NSSolver::~NSSolver(){
    VecDestroy(&b_mom_x);
    VecDestroy(&b_mom_y);
    VecDestroy(&b_pc);
    VecDestroy(&uc);
    VecDestroy(&vc);
    VecDestroy(&pc);
    VecDestroy(&uc_aux);
    VecDestroy(&vc_aux);
    VecDestroy(&pc_aux);
    VecDestroy(&ap);
    
    VecDestroy(&u_face);
    VecDestroy(&v_face);
    VecDestroy(&p_face);
    VecDestroy(&mdotf);
    VecDestroy(&mdotfcorr);
    VecDestroy(&pfcorr);
    VecDestroy(&ucorr);
    VecDestroy(&vcorr);
    VecDestroy(&pcorr);
    VecDestroy(&wf);

    MatDestroy(&A_mom);
    MatDestroy(&A_pc);
    
    KSPDestroy(&ksp);
    
    PetscFinalize();
}

/*
* * Interpolação centrada baseada na distancia dos centroides para o centro da face.
*/
void NSSolver::compute_wf(){
    vector<Cell*> cells = mesh->get_cells();
    vector<Edge*> faces = mesh->get_edges();
    for(PetscInt i = 0; i < faces.size(); i++){
        Edge* face = faces[i];
        int idf = face->id;
        if(!face->is_boundary_face()){
            pair<int,int> nodes_share_face = face->get_link_face_to_cell();
            
            int ic1 = nodes_share_face.first;
            int ic2 = nodes_share_face.second;
            
            Cell* cell1 = cells[ic1];
            Cell* cell2 = cells[ic2];
            
            double d1 = distance(cell1->get_centroid().first, cell1->get_centroid().second, face->get_middle().first, face->get_middle().second);
            double d2 = distance(cell2->get_centroid().first, cell2->get_centroid().second, face->get_middle().first, face->get_middle().second);

            VecSetValue(wf, idf, d2/(d1+d2), INSERT_VALUES);
        }
    }
}

void NSSolver::mom_links_and_sources(double lambda_v){
    vector<Cell*> cells = mesh->get_cells();
   
    MatZeroEntries(A_mom);
    PetscScalar pc_val, ap_val, wf_val, u_val, v_val;
    /*
    * * Convecção-difusão.
    */
    for(PetscInt i = 0; i < cells.size(); ++i){
        // * warm-up
        Cell* c = cells[i];
        int ic = c->id;
        vector<int> &nsigns = c->get_nsigns();
        vector<Edge*> faces = c->get_edges();
        
        // zera tudo
        VecSetValue(b_mom_x, ic, 0.0, INSERT_VALUES);
        VecSetValue(b_mom_y, ic, 0.0, INSERT_VALUES);
        
        for(PetscInt j = 0; j < faces.size(); ++j){
            // * warm-up
            Edge* face = faces[j];
            int idf = face->id;
            int nsign = nsigns[j];
   
            double Df = (mu * face->get_length()) / face->get_df();
            PetscScalar mf;
            VecGetValues(mdotf, 1, &idf, &mf);
            mf *= nsign;           
            
            if(face->is_boundary_face()){
                // * Se uma face é de contorno, então ela tem valor prescrito de velocidade => Mas, a velocidade produto escalar normal nos contornos é sempre zero, então o termo advectivo cancela => Soma difusão no aP e passsa pro outro lado difusão pro termo fonte.
                
                MatSetValue(A_mom, ic, ic, Df, ADD_VALUES);

                VecGetValues(u_face, 1, &idf, &u_val);
                VecGetValues(v_face, 1, &idf, &v_val);

                VecSetValue(b_mom_x, ic, u_val * Df, ADD_VALUES);
                VecSetValue(b_mom_y, ic, v_val * Df, ADD_VALUES);
            }else{
                int N = get_neighbor(face->get_link_face_to_cell(), ic);

                MatSetValue(A_mom, ic, ic, Df + max(mf, 0.0), ADD_VALUES);
                MatSetValue(A_mom, ic, N, -Df - max(-mf,0.0), INSERT_VALUES);
            }
        }

        //* Após o cálculo, eu tenho que aplicar o lambda e salvar o coeficiente da diagonal para uso posterior.
        PetscScalar diag_val;
        MatGetValues(A_mom, 1, &ic, 1, &ic, &diag_val);
        diag_val /= lambda_v;
        MatSetValue(A_mom, ic, ic, diag_val, INSERT_VALUES);
        VecSetValue(ap, ic, diag_val, INSERT_VALUES);
    }

    /*
    * * Calcula pressão nas faces
    */
    vector<Edge*> faces = mesh->get_edges();
    PetscScalar val;
    for(PetscInt i = 0; i < faces.size(); i++){
        Edge* face = faces[i];
        int idf = face->id;
        pair<int,int> id_nodes_share_face = face->get_link_face_to_cell();

        int ic1 = id_nodes_share_face.first;
        int ic2 = id_nodes_share_face.second;

        if(face->is_boundary_face()){
            int neighbor = ic1 != -1 ? ic1 : ic2;
            
            VecGetValues(pc, 1, &neighbor, &val); 
            VecSetValue(p_face, idf, val, INSERT_VALUES);
        }else{
            PetscScalar wf_val, pc1, pc2;
            VecGetValues(wf, 1, &idf, &wf_val);
            VecGetValues(pc, 1, &ic1, &pc1);
            VecGetValues(pc, 1, &ic2, &pc2);
            VecSetValue(p_face, idf, wf_val * pc1 + (1-wf_val) * pc2, INSERT_VALUES);
        }
    }

    /**
     * * Calcula fonte do gradiente de pressão.
     */
    double regularization = ( 1- lambda_v);
    for(PetscInt i = 0; i < cells.size(); ++i){
        Cell* c = cells[i];
        int ic = c->id;
        vector<Edge*> faces = c->get_edges();
        vector<int> &nsigns = c->get_nsigns();
        
        for(PetscInt j = 0; j < faces.size(); ++j){
            int nsign = nsigns[j];
            Edge* face = faces[j];
            int idf = face->id;

            PetscScalar p_face_val;

            VecGetValues(p_face, 1, &idf, &p_face_val);

            VecSetValue(b_mom_x, ic, -p_face_val * face->get_normal().first * nsign * face->get_length(), ADD_VALUES);

            VecSetValue(b_mom_y, ic, -p_face_val * face->get_normal().second * nsign * face->get_length(), ADD_VALUES);
        }

        // & Adiciona contribuição de regularização
        PetscScalar ap_val, uc_val, vc_val;
        VecGetValues(ap, 1, &ic, &ap_val);
        VecGetValues(uc, 1, &ic, &uc_val);
        VecGetValues(vc, 1, &ic, &vc_val);

        VecSetValue(b_mom_x, ic, regularization * ap_val * uc_val, ADD_VALUES);
        VecSetValue(b_mom_y, ic, regularization * ap_val * vc_val, ADD_VALUES);
    }
}

void NSSolver::solve_x_mom(){
    MatAssemblyBegin(A_mom, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A_mom, MAT_FINAL_ASSEMBLY);
    VecAssemblyBegin(b_mom_x);
    VecAssemblyEnd(b_mom_x);

    KSPSetOperators(ksp, A_mom, A_mom);
    KSPSetFromOptions(ksp);
    KSPSolve(ksp, b_mom_x, uc);
}


void NSSolver::solve_y_mom(){
    MatAssemblyBegin(A_mom, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A_mom, MAT_FINAL_ASSEMBLY);
    VecAssemblyBegin(b_mom_y);
    VecAssemblyEnd(b_mom_y);

    KSPSetOperators(ksp, A_mom, A_mom);
    KSPSetFromOptions(ksp);
    KSPSolve(ksp, b_mom_y, vc);
}


// /**
//  * * Realiza a interpolação de momento para obter a velocidade nas faces.
//  */
// void NSSolver::face_velocity(){
//     vector<Edge*> faces = mesh->get_edges();
//     vector<Cell*> cells = mesh->get_cells();
    
//     for(PetscInt i = 0; i < faces.size(); ++i){
//         Edge* face = faces[i];
//         int idf = face->id;

//         if(face->is_boundary_face()) continue; // IGNORE BC
//         else{
//             pair<int,int> nodes_share_face = face->get_link_face_to_cell();
//             int c1 = nodes_share_face.first;
//             int c2 = nodes_share_face.second;

//             double velf_x = wf[idf]*uc[c1] + (1-wf[idf]) * uc[c2];
//             double velf_y = wf[idf]*vc[c1] + (1-wf[idf]) * vc[c2];
//             double vel_f_i = velf_x * face->get_normal().first + velf_y * face->get_normal().second;

//             double v0dp0_x = 0;
//             double v0dp0_y = 0;

//             vector<Edge*> facesOfc1 = cells[c1]->get_edges();
//             vector<int> nsignsc1 = cells[c1]->get_nsigns();
//             for(PetscInt j = 0; j < facesOfc1.size(); j++){
//                 // faces of cell O
//                 Edge* faceC1 = facesOfc1[j];
//                 int nsign = nsignsc1[j];
//                 v0dp0_x = v0dp0_x + p_face[faceC1->id] * faceC1->get_length() * faceC1->get_normal().first * nsign;
//                 v0dp0_y = v0dp0_y + p_face[faceC1->id] * faceC1->get_length() * faceC1->get_normal().second * nsign;
//             }

//             double v1dp1_x = 0;
//             double v1dp1_y = 0;

//             vector<Edge*> facesOfc2 = cells[c2]->get_edges();
//             vector<int> nsignsc2 = cells[c2]->get_nsigns();
//             for(PetscInt j = 0; j < facesOfc2.size(); j++){
//                 // faces of cell N
//                 Edge* faceC2 = facesOfc2[j];
//                 int nsign = nsignsc2[j];
//                 v1dp1_x = v1dp1_x + p_face[faceC2->id] * faceC2->get_length() * faceC2->get_normal().first * nsign;
//                 v1dp1_y = v1dp1_y + p_face[faceC2->id] * faceC2->get_length() * faceC2->get_normal().second * nsign;
//             }
            
//             velf_x = wf[idf] * v0dp0_x/ap[c1] + (1-wf[idf]) * v1dp1_x/ap[c2];
//             velf_y = wf[idf] * v0dp0_y/ap[c1] + (1-wf[idf]) * v1dp1_y/ap[c2];
            
//             double velf_p = velf_x * face->get_normal().first + velf_y * face->get_normal().second;
            
//             double vdotn = vel_f_i + velf_p - (wf[idf]*cells[c1]->get_area()/ap[c1] + (1-wf[idf])*cells[c2]->get_area()/ap[c2]) * (pc[c2] - pc[c1])/face->get_df();
            
//             VecSetValue(mdotf, face->id, vdotn * rho * face->get_length(), INSERT_VALUES);
//         }
//     }
// }

// /**
//  * * Calcula a correção da pressão (p')
//  */
// void NSSolver::solve_pp(){

//     /*Calcula termo fonte*/
//     vector<Cell*> cells = mesh->get_cells();
//     for(int i = 0; i < cells.size(); ++i){
        
//         // * warm-up
//         Cell* c = cells[i];
//         int ic = c->id;
//         vector<Edge*> faces = c->get_edges();
//         vector<int> nsigns = c->get_nsigns();
        
//         VecSetValue(b_pc, ic, 0.0, INSERT_VALUES);
//         for(int j = 0; j < faces.size(); ++j){
//             //* warm-up
//             Edge* face = faces[j];
//             int idf = face->id;
//             int nsign = nsigns[j];

//             VecSetValue(b_pc, ic, -mdotf[idf] * nsign, ADD_VALUES);
//         }
//     }
        
//    // Calcula coeficientes da A
//    MatZeroEntries(A_pc);
//    for(PetscInt i = 0; i < cells.size(); ++i){
//         Cell* c = cells[i];
//         int ic = c->id;

//         vector<Edge*> faces = c->get_edges();
//         vector<int> &nsigns = c->get_nsigns();

//         for(PetscInt j = 0; j < faces.size(); ++j){ // loop over faces of cell
//             int sign = nsigns[j];
//             Edge* face = faces[j];
//             int idf = face->id;

//             if(face->is_boundary_face()) continue;
//             else{
//                 int N = get_neighbor(face->get_link_face_to_cell(), ic);
                
//                 MatSetValues(A_pc, ic, ic, rho * face->get_length() * (wf[idf] * cells[ic]->get_area()/ap[ic] + (1-wf[idf]) * cells[N]->get_area()/ap[N])/face->get_df(), ADD_VALUES);

//                 MatSetValues(A_pc, ic, N, -rho * face->get_length() * (wf[idf] * cells[ic]->get_area()/ap[ic] + (1-wf[idf]) * cells[N]->get_area()/ap[N])/face->get_df(), INSERT_VALUES);
//             }
//         }
//     }

//     // A_pc.row(0).zeros();
//     // A_pc(0,0) = 1;
//     // b_pc[0] = 0;

//     MatAssemblyBegin(A_pc, MAT_FINAL_ASSEMBLY);
//     MatAssemblyEnd(A_pc, MAT_FINAL_ASSEMBLY);
//     VecAssemblyBegin(b_pc);
//     VecAssemblyEnd(b_pc);

//     KSPSetOperators(ksp, A_pc, A_pc);
//     KSPSetFromOptions(ksp);
//     KSPSolve(ksp, b_pc, pcorr);
// }


// void NSSolver::uv_correct() {
//     PetscScalar val;

//     // Interpola pfcorr nas faces
//     vector<Edge*> faces = mesh->get_edges();
//     for(PetscInt i = 0; i < faces.size(); i++){
//         Edge* face = faces[i];
//         int idf = face->id;
//         pair<int,int> id_nodes_share_face = face->get_link_face_to_cell();
//         int ic1 = id_nodes_share_face.first;
//         int ic2 = id_nodes_share_face.second;

//         if(face->is_boundary_face()){
//             int neighbor = ic1 != -1 ? ic1 : ic2;
            
//             VecGetValues(pcorr, 1, &neighbor, &val); 
//             VecSetValue(pfcorr, idf, val, INSERT_VALUES);
//         }else{
//             VecSetValues(pfcorr, idf, wf[idf] * pcorr[ic1] + (1-wf[idf]) * pcorr[ic2], INSERT_VALUES);
//         }
//     }
    
//     // Corrigindo valores das velocidades no centro
//     vector<Cell*> cells = mesh->get_cells();
//     uc_aux = uc; // * salva valor antigo
//     vc_aux = vc; // * salva valor antigo
//     for(PetscInt i = 0; i < cells.size(); i++){
//         Cell* c = cells[i];
//         int ic = c->id;
//         ucorr[ic] = 0;
//         vcorr[ic] = 0;
        
//         vector<Edge*> faces = c->get_edges();
//         vector<int> &nsigns = c->get_nsigns();
//         for(int j = 0; j < faces.size(); j++){
//             Edge* face = faces[j];
//             int nsign = nsigns[j];

//             VecSetValue(ucorr, ic, pfcorr[face->id] * face->get_normal().first * nsign * face->get_length(), ADD_VALUES);
//             VecSetValue(vcorr, ic, pfcorr[face->id] * face->get_normal().second * nsign * face->get_length(), ADD_VALUES);
//         }

//         VecSetValues(ucorr, ic, -ucorr[ic]/ap[ic], INSERT_VALUES);
//         VecSetValues(vcorr, ic, -vcorr[ic]/ap[ic], INSERT_VALUES);
        
//         VecSetValues(uc, ic, ucorr[ic], ADD_VALUES);
//         VecSetValues(vc, ic, vcorr[ic], ADD_VALUES);
//     }

//     // Corrigindo valores das velocidades nas faces
//     for(PetscInt i = 0; i < faces.size(); i++){
//         Edge* face = faces[i];
//         int idf = face->id;
//         if(face->is_boundary_face()) continue;
//         else{
//             pair<int,int> nodes_share_face = face->get_link_face_to_cell();
//             int c1 = nodes_share_face.first;
//             int c2 = nodes_share_face.second;

//             double coeff = wf[idf]*cells[c1]->get_area()/ap[c1] + (1-wf[idf]) * cells[c2]->get_area()/ap[c2];
//             VecSetValue(mdotfcorr, idf, rho*coeff*face->get_length()*(pcorr[c1]-pcorr[c2])/face->get_df(), INSERT_VALUES)

//             VecSetValue(mdotf, idf, mdotfcorr[idf], ADD_VALUES);
//         }
//     }
//     Vec diff;
//     VecDuplicate(uc, &diff);
//     VecWAXPY(diff, -1.0, uc, uc_aux);
//     PetscReal norm_inf;
//     VecNorm(diff, NORM_INFINITY, &norm_inf); 
//     cout << "# Convergência de u: " << norm_inf << endl;

//     ecWAXPY(diff, -1.0, vc, vc_aux);
//     PetscReal norm_inf;
//     VecNorm(diff, NORM_INFINITY, &norm_inf); 
//     cout << "# Convergência de v: " << norm_inf << endl;
//     VecDestroy(&diff);
// };

// void NSSolver::pres_correct(double lambda_p) {
//     vector<Cell*> cells = mesh->get_cells();
    
//     pc_aux = pc; // * salva valores antigos
//     for(PetscInt i = 0; i < cells.size(); i++){
//         int ic = cells[i]->id;
//         VecSetValue(pc, ic, lambda_p * pcorr[ic], ADD_VALUES);
//     }

//     Vec diff;
//     VecDuplicate(uc, &diff);
//     VecWAXPY(diff, -1.0, pc, pc_aux);
//     PetscReal norm_inf;
//     VecNorm(diff, NORM_INFINITY, &norm_inf); 
//     cout << "# Convergência de p: " << arma::norm(pc-pc_aux, "inf") << endl;
//     VecDestroy(&diff);
// };


void NSSolver::export_solution(string filename){

    // ofstream vtk_file(filename);

    // if (!vtk_file.is_open()) {
    //     std::cerr << "ERROR: failed to create file." << filename << "'.\n";
    //     exit(1);
    // }

    // /*Informações padrão*/
    // vtk_file << "# vtk DataFile Version 3.0\n";
    // vtk_file << "Malha 2D com triângulos e temperatura\n";
    // vtk_file << "ASCII\n";
    // vtk_file << "DATASET UNSTRUCTURED_GRID\n\n";

    // vector<Node*> nodes = this->mesh->get_nodes();
    // vtk_file << "POINTS " << nodes.size() << " double" << endl;
    // for(int i = 0; i < nodes.size(); i++){
    //     vtk_file << nodes[i]->x << " " << nodes[i]->y << " " << 0 << endl;
    // }
    // vtk_file << endl;

    // vector<Cell*> cells = mesh->get_cells();
    // vtk_file << "CELLS " << cells.size() << " ";
    // int sum = 0; // vai contar a quantidade de valores a serem lidos depois...
    // for(int i = 0; i < cells.size(); i++)
    // {
    //     if(cells[i]->cellType == 2)
    //         sum += 4;
    //     else if(cells[i]->cellType == 3)
    //         sum += 5;
    // }
    // vtk_file << sum << endl;
    
    // for(int i = 0; i < cells.size(); i++){
    //     vector<Node*>& nodesOfCell = cells[i]->get_nodes();
    //     vtk_file << nodesOfCell.size(); 
    //     for(int j = 0; j < nodesOfCell.size(); j++)
    //         vtk_file << " " << nodesOfCell[j]->id;
    //     vtk_file << endl;
    // }
    // vtk_file << endl;

    // vtk_file << "CELL_TYPES " << cells.size() << endl;
    // for(int i = 0; i < cells.size(); i++){
    //     if(cells[i]->cellType == 2)
    //         vtk_file << "5" << endl;
    //     else if(cells[i]->cellType == 3)
    //         vtk_file << "9" << endl;
    // }
    // vtk_file << endl;
    
    // vtk_file << "CELL_DATA " << cells.size() << endl;
    // vtk_file << "VECTORS velocity double\n";
    // for(int i = 0; i < cells.size(); i++){
    //     double u = uc[i];
    //     double v = vc[i];
    //     // * SALVA U E V.
    //     vtk_file << u << " " << v << " 0.0\n";
    // }


    // vtk_file.close();
}