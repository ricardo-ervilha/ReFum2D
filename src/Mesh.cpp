#include "../include/Mesh.h"
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>

Mesh::Mesh(){
    this->elementTypeToNumNodes[1] = 2; //linha
    this->elementTypeToNumNodes[2] = 3; //triângulo
    this->elementTypeToNumNodes[3] = 4; //quadrilatero
    this->elementTypeToNumNodes[4] = 4; //tetraedro
    this->elementTypeToNumNodes[15] = 1; //ponto

    this->geom_type = -1; // inicia com uma flag.
    this->ncells = 0; //inicializa com 0 para ir contando enquanto for lendo da malha.
}

Mesh::~Mesh(){
}

void Mesh::read_mesh(string filepath){
    ifstream file(filepath);

    string line;

    if(file.is_open()){
        while(getline(file, line)){
            line.erase(remove(line.begin(), line.end(), '\r'), line.end()); /*remove caracteres problemáticos*/
            if (line == "$fileFormat") {
                /*Verificação da versão do arquivo.*/
                getline(file, line);
                istringstream iss(line);
                double version;
                iss >> version;
                if(version == 2.2)
                    continue;
                else
                    cout << "ERROR: Unexpected version of .msh file." << endl;
            } else if(line == "$PhysicalNames"){
                /*Verificação dos grupos físicos*/
                getline(file,line);
                istringstream iss(line);
                iss >> this->totalPhysicalGroups;
                
                for(int i = 0; i < this->totalPhysicalGroups; i++){
                    getline(file, line);
                    istringstream iss(line);
                    int dim, id;
                    string name;
                    iss >> dim >> id >> name;
                    name = name.substr(1, name.size() - 2); // remoção das aspas duplas que vem do arqquivo.
                    
                    PhysicalGroup pg = PhysicalGroup(dim, id, name);
                    this->physicalGroups.emplace(id, pg);

                    if(dim > this->geom_type)
                        this->geom_type = dim; //capturando a geometria do problema inteiro (2D ou 3D)
                }
            } else if(line == "$Nodes"){
                /*Verificação dos nós*/
                getline(file,line);
                istringstream iss(line);
                iss >> this->nnodes;
                
                for(int i = 0; i < this->nnodes; i++){
                    getline(file, line);
                    istringstream iss(line);
                    int id;
                    double x, y, z;
                    iss >> id >> x >> y >> z;
                    Node n = Node(x, y, z);
                    this->nodes.push_back(n);
                }

            }else if(line == "$Elements"){
                /*Verificação dos elementos*/
                getline(file,line);
                istringstream iss(line);
                iss >> this->totalElements;

                for(int i = 0; i < this->totalElements; i++){
                    getline(file, line);
                    istringstream iss(line);
                    int elementId;
                    int elementType;
                    int numTags;
                    iss >> elementId >> elementType >> numTags; 
                    int geometricalEntityTag; // Não tem utilidade.
                    int physicalGroupTag;
                    if(numTags == 2){    
                        iss >> physicalGroupTag >> geometricalEntityTag;
                    }else{
                        cout << "ERROR: numTags != 2." << endl;
                        exit(1);
                    }
                    int numNodes = this->elementTypeToNumNodes[elementType];
                    vector<int> nodeIds;
                    for(int i = 0; i < numNodes; i++){
                        int id;
                        iss >> id;
                        nodeIds.push_back(id-1); //subtrai 1 para ter indexação de 0 e bater com o do vector.
                    }
                    Element e = Element(elementType, nodeIds);
                    this->elements.push_back(e);

                    this->physicalGroups[physicalGroupTag].insert_element_id(elementId-1); //subtrai 1 para ter indexação de 0 e bater com o do vector.

                    if(elementType >= 2 && elementType <= 4){
                        this->ncells++; //incrementa em 1 o número de células.
                        this->cells.push_back(elementId-1);
                    }
                }
            }
        }
        file.close();
    }
    else{
        cout << "ERROR: could not open the file." << endl;
        exit(1);
    }

    this->offset = this->totalElements - this->ncells;
}

void Mesh::pre_processing(vector<BoundaryCondition*> *boundaries){
    /*PRIMEIRA COISA: Computa o conjunto de todas as faces, e para cada célula, já guarda quais faces pertencem a ela.*/
    /*
    Estratégia: -iterar células
                -ordena a face
                -insere em um unordered set
                -unordered set properties...
                    -Insert an element	O(1) (average)
                    -Delete an element	O(1) (average)
                    -Access element by position O(n)
                    -Find element by value	O(1) (average)
                    -Traverse the set	O(n)
    */
    unordered_set<pair<int,int>, CantorPairingFunction> facesUS; //para garantir faces ÚNICAS
    vector<Edge> faces; //vector com as faces.

    unordered_map<pair<int,int>, int, CantorPairingFunction> facesUM; //para recuperar o índice dos elementos já inseridos
    
    int qtdFaces = 0;
    for(int i = 0; i < this->ncells; i++){
        int idCell = this->cells[i]; //converte id do loop para id da célula
        Element& cell = this->elements[idCell]; //recupera a celula como Element

        vector<int>& nodeIds = cell.get_nodes(); //obtém o id dos nodes dessa célula

        /*separa o id dos nós que compõem o polígono. Pela leitura do gmsh, a->b b->c, c->d, ... já estão no sentido anti-horário*/
        vector<Edge> facesToIterate;
        for(int j = 1; j < nodeIds.size(); j++)
            facesToIterate.push_back(Edge(nodeIds[j-1], nodeIds[j]));
        facesToIterate.push_back(Edge(nodeIds[nodeIds.size()-1], nodeIds[0]));

        for (const Edge& face : facesToIterate) {
            pair<int,int> key = face.asKey(); // chave ordenada

            if (facesUS.insert(key).second) {
                facesUM.emplace(key, qtdFaces);
                faces.push_back(face); // mantém ordem original aqui
                cell.insert_face(qtdFaces);
                qtdFaces++;
            } else {
                int idx = facesUM[key];
                cell.insert_face(idx);
            }
        }

   }

    facesUM.clear(); //limpa
    vector<pair<int, int>> facesCorrigido;
    for(int i = 0; i < faces.size(); i++){
       facesCorrigido.push_back(make_pair(faces[i].from, faces[i].to));
       facesUM.emplace(make_pair(faces[i].from, faces[i].to), i);
    }

   this->nfaces = qtdFaces;
   this->faces = facesCorrigido; 
   this->facesPairToKey = facesUM;

   /*Calcula normal, área da face, ponto médio, etc.*/
   for(int i = 0; i < this->nfaces; i++){
        int n1Id = this->faces[i].first;
        int n2Id = this->faces[i].second;
        Node n1 = this->nodes[n1Id];
        Node n2 = this->nodes[n2Id];

        double A_f = hypot(n1.get_x() - n2.get_x(), n1.get_y() - n2.get_y());

        this->faceAreas.push_back(A_f);

        // calcula o ponto médio da face
        double xm = (n1.get_x() + n2.get_x())/2.0;
        double ym = (n1.get_y() + n2.get_y())/2.0;
        double zm = (n1.get_z() + n2.get_z())/2.0;
        Node middle = Node(xm, ym, zm);
        this->faceMiddlePoints.push_back(middle);

        //calculando primeiro os vetores tangentes.
        double tx = (n2.get_x() - n1.get_x())/A_f;
        double ty = (n2.get_y() - n1.get_y())/A_f;

        if(fabs(tx) < 1e-12)
            tx = 0;
        if(fabs(ty) < 1e-12)
            ty = 0;

        // calcula normal
        tuple<double, double> normal = make_tuple(ty, -tx);
        this->normals.push_back(normal);
   }

    /*Cálculo dos centroides & link_face_to_cell*/
    pair<int, int> default_pair(-1,-1);
    vector<pair<int,int>> lftc(this->nfaces, default_pair); //inicializa o vector com tamanho nfaces e default_pair com as flags -1
    this->link_face_to_cell = lftc;
    
    for(int i = 0; i < this->ncells; i++){
        int idCell = this->cells[i];
        Element& cell = this->elements[idCell];
        vector<int>& nodesFromTheCell = cell.get_nodes();

        double xc = 0, yc = 0, zc = 0;
        for(int j = 0; j < nodesFromTheCell.size(); j++){
            Node& n = this->nodes[nodesFromTheCell[j]]; // recupera o nó do polígono 
            xc += n.get_x();
            yc += n.get_y();
            zc += n.get_z();
        }
        xc = xc / nodesFromTheCell.size();
        yc = yc / nodesFromTheCell.size();
        zc = zc / nodesFromTheCell.size();

        Node nc = Node(xc, yc, zc);
        this->centroids.push_back(nc);

        //-------------------------
        //calcula link_face_to_cell
        vector<int>& facesFromCell = cell.get_face_ids();
        for(int j = 0; j < facesFromCell.size(); j++){
            int id = facesFromCell[j];
            if(this->link_face_to_cell[id].first == -1) //verifica se primeiro já foi preenchido
                this->link_face_to_cell[id].first = idCell;
            else //se primeiro já foi preenchido, segundoe estará vazio
                this->link_face_to_cell[id].second = idCell;
        }
    }

    /*Determinação do sinal das normais de cada célula*/
    for(int i = 0; i < this->ncells; i++){ // PARA CADA CÉLULA
        int cellId = this->cells[i]; //conversão do loop index para cell index
        Element& cell = this->elements[cellId]; //recupera célula
        vector<int>& faceIds = cell.get_face_ids(); 
        for(int j = 0; j < faceIds.size(); j++){ //PARA CADA FACE DA CÉLULA
            int faceId = faceIds[j];
            int ic1 = link_face_to_cell[faceId].first; //pega a primeira célula da faceId
            if(cellId == ic1){
                cell.insert_normal_sign(1); // já está apontando OUT
            }
            else if(cellId != -1){
                cell.insert_normal_sign(-1); // estava apontando IN, inverte
            }
            else
            {
                cout << "ERROR: Line 272 in mesh.cpp (problem in normal signs)!" << endl;
                exit(1);
            }
        }
    }

    /*Cálculo do volume das células...*/
    for(int i = 0; i < this->ncells; i++){ // PARA CADA CÉLULA
        int cellId = this->cells[i]; //conversão do loop index para cell index
        Element& cell = this->elements[cellId]; //recupera célula
        double volume = 0;
        vector<int>& faceIds = cell.get_face_ids(); 
        vector<int>& normalSigns = cell.get_normal_sign();
        for(int j = 0; j < faceIds.size(); j++){ //PARA CADA FACE DA CÉLULA
            int gface = faceIds[j];
            //equação 7.39
            volume += get<0>(this->normals[gface]) * normalSigns[j] * this->faceMiddlePoints[gface].get_x() * this->faceAreas[gface];
        }
        this->volumes.push_back(volume);
    }

    /*Calculando o df*/
    int offset = this->totalElements - this->ncells;

    for(int i = 0; i < this->nfaces; i++){ // PARA CADA FACE
        pair<int,int> cells = this->link_face_to_cell[i];
        int cell1 = cells.first - offset; //ids globais 
        int cell2 = cells.second - offset; //ids globais
        
        if(cell1 >= 0 && cell2 >= 0){
            Node* c1 = &this->centroids[cell1]; // centroid da célula 1
            Node* c2 = &this->centroids[cell2]; // centroid da célula 2
            
            //calculando vetor If
            double Ifx = c2->get_x() - c1->get_x();
            double Ify = c2->get_y() - c1->get_y();
            
            tuple<double, double> normal = this->normals[i]; //obtém as componentes da normal daquela face

            double df = Ifx * get<0>(normal) + Ify * get<1>(normal);
            this->deltafs.push_back(df);
        }else{
            this->deltafs.push_back(-1);
        }
    }

    /*Calculando para todos os nós as distâncias deles aos centroides*/
    for(int i = 0; i < this->ncells; i++){ //para cada célula
        int cellId = this->cells[i];
        Element& cell = this->elements[cellId]; //TENHO A CÉLULA
        vector<int>& nodesOfTheCell = cell.get_nodes(); // obtenho ID dos nodes da célula;
        Node& centroidFromCell = this->centroids[i]; //obtém o centroide da célula
        for(int j = 0; j < nodesOfTheCell.size(); j++){ // para cada nó da célula
            int gnode = nodesOfTheCell[j]; //converte em índice global
            Node& n = this->nodes[gnode]; // obtém o nó
            double dist = sqrt(pow(n.get_x() - centroidFromCell.get_x(), 2) + pow(n.get_y() - centroidFromCell.get_y(), 2) + pow(n.get_z() - centroidFromCell.get_z(), 2));
            n.insert_distance_centroids(dist);// armazenar dist no vetor de distcentroid do nó
            n.insert_id_relative_to_centroid(cellId);// armazenar o id da celula lá também
        }
    }

    /*Fazendo a parte das faces de contorno*/
    /*ASSUMINDO QUE VIRÁ NA ORDEM: DOWN, RIGHT, TOP, LEFT e FLUID. (REVER ISSO DEPOIS...)*/
    int contBFaces = 0;
    vector<int> aux_link_face_to_bface(this->nfaces, -1);
    this->link_face_to_bface = aux_link_face_to_bface;
    int contPG = 0;
    for(auto it = this->physicalGroups.begin(); it != this->physicalGroups.end() && it->second.get_dimension() != 2; ++it){ // para cada grupo físico
        vector<int>& elements = it->second.get_element_ids();
        for(int j = 0; j < elements.size(); j++){ // para cada id de elemento
            int idElement = elements[j];
            vector<int>& nodes = this->elements[idElement].get_nodes(); //teoricamenmte é para ser 1D
            int n1 = nodes[0];
            int n2 = nodes[1];
            pair<int, int> faceToTest = make_pair(n1,n2);
            int key = this->facesPairToKey[faceToTest]; //retorna o índice GLOBAL dessa face
            this->link_bface_to_face.push_back(key);
            
            Node* v1 = this->get_node(n1);
            Node* v2 = this->get_node(n2);
            Node* middle = this->get_middle_face(key);
            
            double d1 = sqrt(pow(v1->get_x() - middle->get_x(), 2) + pow(v1->get_y() - middle->get_y(), 2));
            double d2 = sqrt(pow(v2->get_x() - middle->get_x(), 2) + pow(v2->get_y() - middle->get_y(), 2));
            double boundary1 = (*boundaries)[contPG]->apply(v1->get_x(), v1->get_y());
            double boundary2 = (*boundaries)[contPG]->apply(v2->get_x(), v2->get_y());
            double interpolation = (boundary1*d2 + boundary2*d1)/(d1+d2);
            
            this->phib.push_back(interpolation);
            
            this->link_face_to_bface[key] = contBFaces;
            contBFaces++;
        }
        contPG++;
    }

    this->nbfaces = contBFaces;
}

void Mesh::print_deltafs(){
    cout << "#--------------------------------DeltaF's-----------------------------------#" << endl;
    for(int i = 0; i < this->nfaces; i++){ // PARA CADA FACE
        cout << "[ " << i << " ] : " << this->deltafs[i] << endl;
    }
}

void Mesh::print_info_geral(){
    cout << "#------------------------------Infos. Gerais-----------------------------------------#" << endl;
    cout << "Geometria do problema: " << this->geom_type << "D" << endl;
    cout << endl;
}

void Mesh::print_info_nodes(){
    cout << "#----------------------------------------Nós-----------------------------------------#" << endl;
    cout << "Quantidade total de nós: " << this->nnodes << endl;
    cout << "Lista de coordenadas:" << endl;
    for(int i = 0; i < this->nnodes; i++)
        cout << "[" << i << "]" << "\tx: " << this->nodes[i].get_x() << "\ty:" << this->nodes[i].get_y() << "\tz:" << this->nodes[i].get_z() << endl;
    cout << endl;
}

void Mesh::print_info_elements(){
    cout << "#----------------------------------------Elementos-----------------------------------#" << endl;
    cout << "Quantidade de elementos: " << this->totalElements << endl;
    cout << "Informações dos elementos:" << endl;
    for(int i = 0; i < this->totalElements; i++){
        cout << "[" << i << "]" << "\tTipo: " << this->elements[i].get_element_type() << " \tNós que o compõem: ";
        vector<int>& ns = this->elements[i].get_nodes();
        for(int j = 0; j < ns.size(); j++)
            cout << ns[j] << "   ";
        cout << endl;
    } 
    cout << endl;
}

void Mesh::print_info_physical_groups(){
    cout << "#-----------------------------Grupos Físicos----------------------------------------#" << endl;
    cout << "Quantidade de grupos físicos: " << this->totalPhysicalGroups << endl;
    cout << "Informação dos grupos físicos:" << endl;
    for(auto it = this->physicalGroups.begin(); it != this->physicalGroups.end(); ++it){
        cout << "[" << it->second.get_id() << "]" <<  "\tDimensão: " << it->second.get_dimension() << " Nome: " << it->second.get_name() << "   Id dos elementos que pertencem a ele: ";
        vector<int>& es = it->second.get_element_ids();
        for(int i = 0; i < es.size(); i++)
        cout <<  es[i] << " ";
        cout << endl;
    }
    cout << endl;
}

void Mesh::print_info_cells(){
    cout << "#------------------------------------Células----------------------------------------#" << endl;
    cout << "Quantidade de células: " << this->ncells << endl;
    cout << "Id dos elementos que são células: ";
    for(int i = 0; i < this->ncells; i++)
        cout << this->cells[i] << " ";
    cout << endl;
    cout << endl;
}

void Mesh::print_info_faces(){
    cout << "#------------------------------------Faces------------------------------------------#" << endl;
    cout << "Quantidade de faces: " << this->nfaces << endl;
    cout << "Id das faces e seus nós constituintes:" << endl;
    for(int i = 0; i < this->nfaces; i++)
        cout  << "[" << i << "] " << this->faces[i].first << " " << this->faces[i].second << endl;
    cout << "#-------------------------------------Faces das Células-----------------------------#" << endl;
    
    cout << endl;
    cout << "Células e suas faces: " << endl;
    for(int i = 0; i < this->ncells; i++){
        int cellId = this->cells[i]; //casa a indexação do loop com a indexação verdadeira das células
        vector<int>& faceIds = this->elements[cellId].get_face_ids();
        cout << "[" << cellId << "] ";
        for(int j = 0; j < faceIds.size(); j++)
            cout <<  faceIds[j] << " ";
        cout << endl;
    }
    cout << endl;

    cout<< "#------------------------------------Áreas das Faces-----------------------------------------#" << endl;
    for(int i = 0; i < this->nfaces; i++)
        cout  << "[" << i << "] " << this->faceAreas[i] << endl;
    cout << endl;
}

void Mesh::print_info_centroids(){  
    cout << "#-------------------------------------Centroides-------------------------------------#" << endl;
    for(int i = 0; i < this->ncells; i++){
        int cellId = this->cells[i]; //casa a indexação do loop com a indexação verdadeira das células
        cout << "[" << cellId << "] \tx:" << this->centroids[i].get_x() << " \ty: " << this->centroids[i].get_y() << " \tz: " << this->centroids[i].get_z() << endl;
    }
    cout << endl;
}

void Mesh::print_info_link_face_to_cell(){
    cout << "#-------------------------------------Link Face to Cell------------------------------#" << endl;
    for(int i = 0; i < this->nfaces; i++){
        cout <<  "[" << i << "] " << this->link_face_to_cell[i].first << " " << this->link_face_to_cell[i].second << endl;
    }
    cout << endl;
}

void Mesh::print_info_normal_signs(){
    cout << "#----------------------------------------Sinal das Normais---------------------------#" << endl;
    for(int i = 0 ; i < this->ncells; i++){
        int cellId = this->cells[i]; 
        Element& cell = this->elements[cellId];
        vector<int>& normalSigns = cell.get_normal_sign();
        vector<int>& faceIds = cell.get_face_ids();
        for(int j=0; j < normalSigns.size(); j++)
            cout << "[" << faceIds[j] << "] -> " << normalSigns[j] << " ";
        cout << endl;
    }
    cout << endl;
}

void Mesh::print_info_normal_values(){
    cout << "#--------------------------------------------Normais----------------------------------#" << endl;
    for(int i = 0; i < this->nfaces; i++){
        cout << "[" << i << "]" << "\tx: " << get<0>(normals[i]) << "\ty: " << get<1>(normals[i]) << endl;
    }
    cout << endl;
}

void Mesh::print_info_volumes(){
    cout << "#--------------------------------------------Volumes----------------------------------#" << endl;
    for(int i = 0; i < this->ncells; i++){
        int cellId = this->cells[i]; 
        cout << "Volume da celula " << cellId << " : " << volumes[i] << endl;
    }
    cout << endl;
}

void Mesh::print_info_link_face_to_bface(){
    cout << "#------------------------------------Link Face to Bface-------------------------------#" << endl;
    for(int i = 0; i < this->link_face_to_bface.size(); i++){
        cout << "[" << i << "] = " << this->link_face_to_bface[i] << endl;
    }    
    cout << endl;
}

void Mesh::print_info_link_bface_to_face(){  
    cout << "#--------------------------------Link Bface to Face-----------------------------------#" << endl;
    for(int i = 0 ; i < this->link_bface_to_face.size(); i++){
        cout << "[" << i << "] = " << this->link_bface_to_face[i] << endl;
    }
    cout << endl;
}

void Mesh::print_info_distance_node_to_centroids(){
    cout << "#-----------------------------------Distancia aos Centroides -------------------------#" << endl;
    for(int i = 0; i < this->nnodes; i++){
        vector<double>& distCentroides = this->nodes[i].get_distance_centroids();
        vector<int>& idCellsOfCentroids = this->nodes[i].get_id_relative_to_centroids();
        for(int j = 0; j < distCentroides.size(); j++)
            cout << "Node: " << i << " \tCélula: " << idCellsOfCentroids[j] << " \tDistance: " << distCentroides[j] << endl;
    }
    cout << endl;
}

void Mesh::mesh_summary(){
    print_info_geral();
    print_info_nodes();
    print_info_elements();
    print_info_physical_groups();
    print_info_cells();
    print_info_faces();
    print_info_centroids();
    print_info_link_face_to_cell();
    print_info_normal_signs();
    print_info_normal_values();
    print_info_volumes();
    print_info_link_face_to_bface();
    print_info_link_bface_to_face();
    print_deltafs();
    print_info_distance_node_to_centroids();
}