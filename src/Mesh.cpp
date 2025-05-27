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
    vector<pair<int,int>> faces; //vector com as faces.

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

        /*Criando os pairs ordenados usando minmax.*/
        pair<int,int> ab = minmax(a,b); 
        pair<int,int> bc = minmax(b,c);
        pair<int,int> ca = minmax(c,a); 

        vector<pair<int,int>> facesToIterate;
        facesToIterate.push_back(ab);
        facesToIterate.push_back(bc);
        facesToIterate.push_back(ca);

        // iterando para ab, bc, ca
        for(int j = 0; j < facesToIterate.size(); j++){
            if(facesUS.insert(facesToIterate[j]).second){
                // inserção da face no Unordered set deu certo
                facesUM.emplace(facesToIterate[j], qtdFaces); // insere no unordered map
                faces.push_back(facesToIterate[j]); //insere no vector
                cell.insertFace(qtdFaces); // insere a face na célula
                qtdFaces++; // aumenta 1 para dizer que número de faces únicas aumentou
            }else{
                // inserção da face deu errado (já está lá)
                int idx = facesUM[facesToIterate[j]]; // busca O(1)
                cell.insertFace(idx); // insere a face na célula
            }
        }

   }

   this->nfaces = qtdFaces;
   this->faces = faces; 
}

void Mesh::meshSummary(){
    /**/
    cout << "Geometria do problema: " << this->geom_type << "D" << endl;
    cout << "Quantidade total de nós: " << this->nnodes << endl;
    cout << "Lista de coordenadas dos nodes (indexados de 0 até N-1):" << endl;
    for(int i = 0; i < this->nnodes; i++)
        cout << "x: " << this->nodes[i].getX() << "\ty:" << this->nodes[i].getY() << "\tz:" << this->nodes[i].getZ() << endl;
    cout << "#----------------------------------------------------------------------------------------#" << endl;
    cout << "Quantidade de elementos: " << this->totalElements << endl;
    cout << "Informação dos elementos: (indexados de 0 até N-1)" << endl;
    for(int i = 0; i < this->totalElements; i++){
        cout << "Tipo: " << this->elements[i].getElementType() << " | Nós que o compõem: ";
        vector<int>* ns = this->elements[i].getNodes();
        for(int j = 0; j < (*ns).size(); j++)
            cout << (*ns)[j] << " ";
        cout << endl;
    } 
    cout << "#----------------------------------------------------------------------------------------#" << endl;
    cout << "Quantidade de grupos físicos: " << this->totalPhysicalGroups << endl;
    cout << "Informação dos grupos físicos: (indexados pelo que vem do arquivo)" << endl;
    for(auto it = this->physicalGroups.begin(); it != this->physicalGroups.end(); ++it){
        cout << "ID: " << it->second.getId() <<  " Dimensão: " << it->second.getDimension() << " Nome: " << it->second.getName() << " \tId dos elementos que pertencem a ele: ";
        vector<int>* es = it->second.getElementIds();
        for(int i = 0; i < (*es).size(); i++)
        cout <<  (*es)[i] << " ";
        cout << endl;
    }
    cout << "#----------------------------------------------------------------------------------------#" << endl;
    cout << "Quantidade de células: " << this->ncells << endl;
    cout << "Id dos elementos que são células: ";
    for(int i = 0; i < this->ncells; i++)
        cout << this->cells[i] << " ";
    cout << endl;
    cout << "#----------------------------------------------------------------------------------------#" << endl;
    cout << "Quantidade de faces: " << this->nfaces << endl;
    cout << "Id das faces e seus nós constituintes (Indexados de 0 até NumFaces - 1):" << endl;
    for(int i = 0; i < this->nfaces; i++)
        cout  << "[" << i << "] " << this->faces[i].first << " " << this->faces[i].second << endl;
    cout << "#----------------------------------------------------------------------------------------#" << endl;
    cout << "Células e suas faces... (Tudo no anti-horário):" << endl;
    for(int i = 0; i < this->ncells; i++){
        int cellId = this->cells[i]; //casa a indexação do loop com a indexação verdadeira das células
        vector<int>* faceIds = this->elements[cellId].getFaceIds();
        cout << "[" << cellId << "] ";
        for(int j = 0; j < (*faceIds).size(); j++)
            cout <<  (*faceIds)[j] << " ";
        cout << endl;
    }
}