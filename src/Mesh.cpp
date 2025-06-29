#include "../include/Mesh.h"

Mesh::Mesh(){
    this->elementTypeToNumNodes[1] = 2; //Linha tem dois nós
    this->elementTypeToNumNodes[2] = 3; //Triângulo tem três nós
    this->elementTypeToNumNodes[3] = 4; //Quadrilatero tem quatro nós
    this->elementTypeToNumNodes[15] = 1; //Ponto tem um nó

    xmin = numeric_limits<double>::max();
    ymin = numeric_limits<double>::max();

    xmax = numeric_limits<double>::lowest();
    ymax = numeric_limits<double>::lowest();
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
                /*Verificação de versão do arquivo se é 2 ASCII.*/
                getline(file, line);
                istringstream iss(line);
                double version;
                iss >> version;
                if(version == 2.2)
                    continue;
                else
                    cout << "ERROR: Unexpected version of .msh file." << endl;
            }else if(line == "$PhysicalNames"){
                /*Inicia leitura dos grupos físicos*/
                getline(file,line);
                istringstream iss(line);
                int qtdPhysicalGroups;
                iss >> qtdPhysicalGroups;
                
                for(int i = 0; i < qtdPhysicalGroups; i++){
                    getline(file, line);
                    istringstream iss(line);
                    int dim, id;
                    string name;
                    iss >> dim >> id >> name;
                    name = name.substr(1, name.size() - 2); // Remoção das aspas duplas que vem junto com o nome.
                    
                    if(dim == 2){    
                        PhysicalEntity* pg = new PhysicalEntity(id, name);
                        this->physicalentities.emplace(id, pg);
                    }
                }
            } else if(line == "$Nodes"){
                /*Inicia leitura dos nós*/
                getline(file,line);
                istringstream iss(line);
                int nnodes;
                iss >> nnodes;
                
                for(int i = 0; i < nnodes; i++){
                    getline(file, line);
                    istringstream iss(line);
                    int id;
                    double x, y, z;
                    iss >> id >> x >> y >> z;
                    Node* n = new Node(i, x, y);
                    this->nodes.push_back(n);

                    if (x < xmin) xmin = x;
                    if (x > xmax) xmax = x;
                    if (y < ymin) ymin = y;
                    if (y > ymax) ymax = y;
                }
            } else if(line == "$Elements"){
                int idCell = 0; // variável auxiliar para ir pegando o id das células.

                /*Inicia leitura dos elementos*/
                getline(file,line);
                istringstream iss(line);
                int totalElements;
                iss >> totalElements;

                for(int i = 0; i < totalElements; i++){
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
                    vector<Node*> elNodes;
                    
                    for(int j = 0; j < numNodes; j++){
                        int id;
                        iss >> id;
                        elNodes.push_back(this->nodes[id - 1]);
                    }

                    if(elementType == 2 || elementType == 3){
                        // é uma célula.
                        Cell* c = new Cell(idCell, elementType, elNodes);
                        this->cells.push_back(c);
                        idCell++; // incrementa 1
                        this->physicalentities[physicalGroupTag]->add_element(c);
                    }
                }
            }
        }
    }

    /*Após a leitura, chama automaticamente o pre_processing.*/
    this->pre_processing();
}

void Mesh::pre_processing(){
    
    // Quando dizer face, considere a mesma coisa que aresta (Edge).
    /*
        1ª etapa: Computa o conjunto de todas as faces, e para cada célula, já guarda quais faces pertencem a ela.
    */
    /*
    Estratégia para isso:
            -iterar células
                -ordena a face da célula
                -insere em um unordered set
    */
    unordered_set<pair<int,int>, CantorPairingFunction> uniqueFaces;
    vector<pair<int,int>> edgesAux; // guarda as faces (é um auxiliar ao uniqueFaces pois ele não guarda ordem, e preciso disso).
    for(int i = 0; i < this->cells.size(); i++){

        vector<pair<int,int>> FacesOfCell; // Faces dessa célula
        vector<Node*>& nodesOfCell = this->cells[i]->get_nodes(); // Nós dessa céluça
        for(int j = 1; j < nodesOfCell.size(); j++)
            FacesOfCell.push_back(make_pair(nodesOfCell[j-1]->id, nodesOfCell[j]->id));
        FacesOfCell.push_back(make_pair(nodesOfCell[nodesOfCell.size() - 1]->id, nodesOfCell[0]->id));

        // a partir daqui já tenho uma lista com as faces daquela célula.
        // então para cada face dessa célula...
        for(int j = 0; j < FacesOfCell.size(); j++){
            // faço um trick: vou colocar as faces ordenadas.
            int id1 = FacesOfCell[j].first;
            int id2 = FacesOfCell[j].second;
            pair<int,int> ordenered_edge = make_pair(min(id1, id2), max(id1,id2));
            
            // O(1)
            if(uniqueFaces.find(ordenered_edge) != uniqueFaces.end()){
                // está presente em uniqueFaces.
                
                int k;
                for(k = 0; k < edgesAux.size(); k++)
                    if((edgesAux[k].first == id1 && edgesAux[k].second == id2) || (edgesAux[k].second == id1 && edgesAux[k].first == id2))
                        break;
                
                // k agora armazena o índice dessa face.
                cells[i]->insert_edge(this->edges[k]);
            }else {
                // não está presente em uniqueFaces.
                
                // insiro a face na célula na ordem que veio antes. Uso o próprio size para facilitar esse controle.
                Edge* e = new Edge(edgesAux.size(), this->nodes[id1], this->nodes[id2]);
                this->edges.push_back(e); //registro isso também no vetor global de edges. 
                cells[i]->insert_edge(e);
                
                // adiciono em uniquefaces e faces.
                uniqueFaces.insert(ordenered_edge); // insiro essa nova face, registrando uma nova face única.
                edgesAux.push_back(FacesOfCell[j]);  // guardo isso no faces também, pois ele anda junto com o uniqueFaces.
            }
        }
    }

    /*2ª etapa: calcular o link_face_to_cell*/
    for(int i = 0; i < this->cells.size(); i++){
        vector<Edge*>& facesFromCell = cells[i]->get_edges();
        for(int j = 0; j < facesFromCell.size(); j++){
            
            if(facesFromCell[j]->get_link_face_to_cell().first == -1) // verifica pra ver se a posição está vazia.
                facesFromCell[j]->set_link_face_to_cell(cells[i]->id, 1);
            else
                facesFromCell[j]->set_link_face_to_cell(cells[i]->id, 2);
        }
    }

    /*3ª etapa: calcular o sinal das normais*/
    for(int i = 0; i < this->cells.size(); i++){
        vector<Edge*>& facesFromCell = cells[i]->get_edges();
        for(int j = 0; j < facesFromCell.size(); j++){
            int idCell1 = facesFromCell[j]->get_link_face_to_cell().first; // pega a primeira célula do link_face_to_cell
            
            if(idCell1 == cells[i]->id) cells[i]->insert_nsign(1); // já aponta para fora
            else if(idCell1 != -1)  cells[i]->insert_nsign(-1); //conserta, fazendo apontar para fora
            else{
                cout << "ERROR: Line 272 in mesh.cpp (problem in normal signs)!" << endl;
                exit(1);
            }
        }
    }

    /*4º etapa: usar o que foi feito anteriormente para calcular o volume*/
    for(int i = 0; i < this->cells.size(); i++){
        vector<Edge*>& facesFromCell = cells[i]->get_edges();
        double area = 0;
        vector<int>& nsigns = cells[i]->get_nsigns();
        for(int j = 0; j < facesFromCell.size(); j++){
            pair<double,double>& normal = facesFromCell[j]->get_normal();
            pair<double,double>& middle = facesFromCell[j]->get_middle();
            area += normal.first * nsigns[j] * middle.first * facesFromCell[j]->get_length();
        }
        cells[i]->set_area(area);
    }

    /*5ª etapa: calculo do df das arestas.*/
    for(int i = 0; i < this->edges.size(); i++){
        pair<int,int> cellIds = this->edges[i]->get_link_face_to_cell();
        if(this->edges[i]->is_boundary_face()){
            /*caso seja boundary face: calculo da unica celula adjacente e o centro da face. (dB)*/
            int cellAdj = cellIds.first != -1 ? cellIds.first : cellIds.second;
            pair<double,double>& cell = cells[cellAdj]->get_centroid();
            pair<double,double>& mid = this->edges[i]->get_middle();
            double db = distance(cell.first, cell.second, mid.first, mid.second);
            this->edges[i]->set_df(db);
        }else{
            /*caso não seja boundary face: cálculo entre o centroid das 2 celulas adjacentes. (df)*/
            pair<double,double>& c1 = cells[cellIds.first]->get_centroid();
            pair<double,double>& c2 = cells[cellIds.second]->get_centroid();
            double df = distance(c1.first, c1.second, c2.first, c2.second);
            this->edges[i]->set_df(df);
        }
    }
}