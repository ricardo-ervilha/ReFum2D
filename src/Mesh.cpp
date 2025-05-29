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

void Mesh::readMesh(string filepath){
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

                    this->physicalGroups[physicalGroupTag].insertElementId(elementId-1); //subtrai 1 para ter indexação de 0 e bater com o do vector.

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

    this->preProcessing();
}

/*Função hash para o unordered set*/
struct CantorPairingFunction {
    // injetora para inteiros >= 0.
    std::size_t operator()(const std::pair<int, int>& p) const noexcept {
        int a = p.first;
        int b = p.second;
        return (a + b) * (a + b + 1) / 2 + b;
    }
};

struct Edge {
    int from, to; // mantém a ordem original

    Edge(int a, int b) : from(a), to(b) {}

    pair<int, int> asKey() const {
        return (from < to) ? make_pair(from, to) : make_pair(to, from);
    }
};

void Mesh::preProcessing(){
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

        vector<int>* nodeIds = cell.getNodes(); //obtém o id dos nodes dessa célula

        /*separa o id dos três nós que compõem o triângulo. Pela leitura do gmsh, a->b b->c e c->a já estão no sentido anti-horário*/
        int a =  (*nodeIds)[0];
        int b =  (*nodeIds)[1];
        int c =  (*nodeIds)[2];

        Edge ab(a, b);
        Edge bc(b, c);
        Edge ca(c, a);

        vector<Edge> facesToIterate = {ab, bc, ca};

        for (const Edge& face : facesToIterate) {
            pair<int,int> key = face.asKey(); // chave ordenada

            if (facesUS.insert(key).second) {
                facesUM.emplace(key, qtdFaces);
                faces.push_back(face); // mantém ordem original aqui
                cell.insertFace(qtdFaces);
                qtdFaces++;
            } else {
                int idx = facesUM[key];
                cell.insertFace(idx);
            }
        }

   }

    vector<pair<int, int>> facesCorrigido;
    for(int i = 0; i < faces.size(); i++){
        facesCorrigido.push_back(make_pair(faces[i].from, faces[i].to));
    }

   this->nfaces = qtdFaces;
   this->faces = facesCorrigido; 

   /*Calcula normal, área da face, ponto médio, etc.*/
   for(int i = 0; i < this->nfaces; i++){
        int n1Id = this->faces[i].first;
        int n2Id = this->faces[i].second;
        Node n1 = this->nodes[n1Id];
        Node n2 = this->nodes[n2Id];

        // calcula área da face
        double A_f = sqrt(pow(n1.getX() - n2.getX(), 2) + pow(n1.getY() - n2.getY(), 2) + pow(n1.getZ() - n2.getZ(), 2)); 

        this->faceAreas.push_back(A_f);

        // calcula o ponto médio da face
        double xm = (n1.getX() + n2.getX())/2.0;
        double ym = (n1.getY() + n2.getY())/2.0;
        double zm = (n1.getZ() + n2.getZ())/2.0;
        Node middle = Node(xm, ym, zm);
        this->faceMiddlePoints.push_back(middle);

        //calculando primeiro os vetores tangentes.
        double tx = (n2.getX() - n1.getX())/A_f;
        double ty = (n2.getY() - n1.getY())/A_f;

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
        vector<int>* nodesFromTheCell = cell.getNodes();

        Node& n1 = this->nodes[(*nodesFromTheCell)[0]]; 
        Node& n2 = this->nodes[(*nodesFromTheCell)[1]]; 
        Node& n3 = this->nodes[(*nodesFromTheCell)[2]]; 

        double xc = (n1.getX() + n2.getX() + n3.getX())/3.0;
        double yc = (n1.getY() + n2.getY() + n3.getY())/3.0;
        double zc = (n1.getZ() + n2.getZ() + n3.getZ())/3.0;

        Node nc = Node(xc, yc, zc);
        this->centroids.push_back(nc);

        //-------------------------
        //calcula link_face_to_cell
        vector<int>* facesFromCell = cell.getFaceIds();
        for(int j = 0; j < (*facesFromCell).size(); j++){
            int id = (*facesFromCell)[j];
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
        vector<int>* faceIds = cell.getFaceIds(); 
        for(int j = 0; j < (*faceIds).size(); j++){ //PARA CADA FACE DA CÉLULA
            int faceId = (*faceIds)[j];
            int ic1 = link_face_to_cell[faceId].first; //pega a primeira célula da faceId
            if(cellId == ic1){
                cell.insertNormalSign(1); // já está apontando OUT
            }
            else if(cellId != -1){
                cell.insertNormalSign(-1); // estava apontando IN, inverte
            }
            else
            {
                cout << "ERROR: Line 272 in mesh.cpp (problem in normal signs)!" << endl;
                exit(1);
            }
        }
    }

    for(int i = 0; i < this->ncells; i++){ // PARA CADA CÉLULA
        int cellId = this->cells[i]; //conversão do loop index para cell index
        Element& cell = this->elements[cellId]; //recupera célula
        double volume = 0;
        vector<int>* faceIds = cell.getFaceIds(); 
        vector<int>* normalSigns = cell.getNormalSigns();
        for(int j = 0; j < (*faceIds).size(); j++){ //PARA CADA FACE DA CÉLULA
            int gface = (*faceIds)[j];
            //equação 7.39
            volume += get<0>(this->normals[gface]) * (*normalSigns)[j] * this->faceMiddlePoints[gface].getX() * this->faceAreas[gface];
        }
        this->volumes.push_back(volume);
    }

    /*Calculando para todos os nós as distâncias deles aos centroides*/
    for(int i = 0; i < this->ncells; i++){ //para cada célula
        int cellId = this->cells[i];
        Element& cell = this->elements[cellId]; //TENHO A CÉLULA
        vector<int>* nodesOfTheCell = cell.getNodes(); // obtenho ID dos nodes da célula;
        Node& centroidFromCell = this->centroids[i]; //obtém o centroide da célula
        for(int j = 0; j < (*nodesOfTheCell).size(); j++){ // para cada nó da célula
            int gnode = (*nodesOfTheCell)[j]; //converte em índice global
            Node& n = this->nodes[gnode]; // obtém o nó
            double dist = sqrt(pow(n.getX() - centroidFromCell.getX(), 2) + pow(n.getY() - centroidFromCell.getY(), 2) + pow(n.getZ() - centroidFromCell.getZ(), 2));
            n.insertDistanceCentroids(dist);// armazenar dist no vetor de distcentroid do nó
            n.insertIdCellRelativeToCentroid(cellId);// armazenar o id da celula lá também
        }
    }
}

void Mesh::meshSummary(){
    /**/
    cout << "#------------------------Infos. Gerais & Nós-----------------------------------------#" << endl;
    cout << "Geometria do problema: " << this->geom_type << "D" << endl;
    cout << "Quantidade total de nós: " << this->nnodes << endl;
    cout << "Lista de coordenadas dos nodes (indexados de 0 até N-1):" << endl;
    for(int i = 0; i < this->nnodes; i++)
        cout << "x: " << this->nodes[i].getX() << "\ty:" << this->nodes[i].getY() << "\tz:" << this->nodes[i].getZ() << endl;
    cout << "#---------------------------Elementos------------------------------------------#" << endl;
    cout << "Quantidade de elementos: " << this->totalElements << endl;
    cout << "Informação dos elementos: (indexados de 0 até N-1)" << endl;
    for(int i = 0; i < this->totalElements; i++){
        cout << "Tipo: " << this->elements[i].getElementType() << " | Nós que o compõem: ";
        vector<int>* ns = this->elements[i].getNodes();
        for(int j = 0; j < (*ns).size(); j++)
            cout << (*ns)[j] << " ";
        cout << endl;
    } 
    cout << "#-----------------------------Grupos Físicos----------------------------------------------#" << endl;
    cout << "Quantidade de grupos físicos: " << this->totalPhysicalGroups << endl;
    cout << "Informação dos grupos físicos: (indexados pelo que vem do arquivo)" << endl;
    for(auto it = this->physicalGroups.begin(); it != this->physicalGroups.end(); ++it){
        cout << "ID: " << it->second.getId() <<  " Dimensão: " << it->second.getDimension() << " Nome: " << it->second.getName() << " \tId dos elementos que pertencem a ele: ";
        vector<int>* es = it->second.getElementIds();
        for(int i = 0; i < (*es).size(); i++)
        cout <<  (*es)[i] << " ";
        cout << endl;
    }
    cout << "#------------------------------------Células------------------------------------------#" << endl;
    cout << "Quantidade de células: " << this->ncells << endl;
    cout << "Id dos elementos que são células: ";
    for(int i = 0; i < this->ncells; i++)
        cout << this->cells[i] << " ";
    cout << endl;
    cout << "#------------------------------------Faces----------------------------------------------#" << endl;
    cout << "Quantidade de faces: " << this->nfaces << endl;
    cout << "Id das faces e seus nós constituintes (Indexados de 0 até NumFaces - 1):" << endl;
    for(int i = 0; i < this->nfaces; i++)
        cout  << "[" << i << "] " << this->faces[i].first << " " << this->faces[i].second << endl;
    cout << "#-------------------------------------Faces das Células----------------------------------------#" << endl;
    cout << "Células e suas faces... (Tudo no anti-horário):" << endl;
    for(int i = 0; i < this->ncells; i++){
        int cellId = this->cells[i]; //casa a indexação do loop com a indexação verdadeira das células
        vector<int>* faceIds = this->elements[cellId].getFaceIds();
        cout << "[" << cellId << "] ";
        for(int j = 0; j < (*faceIds).size(); j++)
            cout <<  (*faceIds)[j] << " ";
        cout << endl;
    }
    cout << "#-------------------------------------Centroides----------------------------------------#" << endl;
    for(int i = 0; i < this->ncells; i++){
        int cellId = this->cells[i]; //casa a indexação do loop com a indexação verdadeira das células
        cout << "[" << cellId << "] \tx:" << this->centroids[i].getX() << " \ty: " << this->centroids[i].getY() << " \tz: " << this->centroids[i].getZ() << endl;
    }
    cout << "#-------------------------------------Link face to cell----------------------------------------#" << endl;
    for(int i = 0; i < this->nfaces; i++){
        cout <<  "[" << i << "] " << this->link_face_to_cell[i].first << " " << this->link_face_to_cell[i].second << endl;
    }
    cout << "#----------------------------------------Sinal das Normais------------------------------------------#" << endl;
    for(int i = 0 ; i < this->ncells; i++){
        int cellId = this->cells[i]; 
        Element& cell = this->elements[cellId];
        vector<int>* normalSigns = cell.getNormalSigns();
        for(int j=0; j < (*normalSigns).size(); j++)
            cout << (*normalSigns)[j] << " ";
        cout << endl;
    }
    cout << "#--------------------------------------------Normais------------------------------------------------#" << endl;
    for(int i = 0; i < this->nfaces; i++){
        cout <<  "x: " << get<0>(normals[i]) << " \ty: " << get<1>(normals[i]) << endl;
    }
    cout << "#--------------------------------------------Volumes------------------------------------------------#" << endl;
    for(int i = 0; i < this->ncells; i++){
        int cellId = this->cells[i]; 
        cout << "Volume da celula " << cellId << " : " << volumes[i] << endl;
    }
    cout << "#--------------------------------------------Distancia aos Centroides ------------------------------------------------#" << endl;
    for(int i = 0; i < this->nnodes; i++){
        vector<double>* distCentroides = this->nodes[i].getDistanceCentroids();
        vector<int>* idCellsOfCentroids = this->nodes[i].getIdCellRelativeToCentroid();
        for(int j = 0; j < (*distCentroides).size(); j++)
            cout << "Node: " << i << " \tCélula: " << (*idCellsOfCentroids)[j] << " \tDistance: " << (*distCentroides)[j] << endl;
    }
}