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
                getline(file, line);
                istringstream iss(line);
                double version;
                iss >> version;
                if(version == 2.2)
                    continue;
                else
                    cout << "ERROR: Unexpected version of .msh file." << endl;
            } else if(line == "$PhysicalNames"){
                getline(file,line);
                istringstream iss(line);
                iss >> this->totalPhysicalEntities;
                
                for(int i = 0; i < this->totalPhysicalEntities; i++){
                    getline(file, line);
                    istringstream iss(line);
                    int dim, id;
                    string name;
                    iss >> dim >> id >> name;
                    name = name.substr(1, name.size() - 2);
                    PhysicalGroup pg = PhysicalGroup(dim, id, name);
                    this->physicalGroups.emplace(id, pg);

                    if(dim > this->geom_type)
                        this->geom_type = dim; //capturando a geometria do problema inteiro 
                }
            } else if(line == "$Nodes"){
                getline(file,line);
                istringstream iss(line);
                iss >> this->nnodes;
                
                for(int i = 0; i < this->nnodes; i++){
                    getline(file, line);
                    istringstream iss(line);
                    int id;
                    double x, y, z;
                    iss >> id >> x >> y >> z;
                    Node n = Node(id, x, y, z);
                    this->nodes.emplace(id, n);
                }
            }else if(line == "$Elements"){
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
                    int geometricalEntityTag; // ainda não tenho certeza onde usar isso.
                    int physicalGroupTag;
                    if(numTags == 2){    
                        iss >> physicalGroupTag >> geometricalEntityTag;
                    }
                    int numNodes = this->elementTypeToNumNodes[elementType];
                    vector<int> nodeIds;
                    for(int i = 0; i < numNodes; i++){
                        int id;
                        iss >> id;
                        nodeIds.push_back(id);
                    }
                    Element e = Element(elementId, elementType, nodeIds);
                    this->elements.emplace(elementId, e);

                    this->physicalGroups[physicalGroupTag].insertElementId(elementId);

                    if(elementType >= 2 && elementType <= 4)
                        this->ncells++; //incrementa em um o número de células.
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

void Mesh::preProcessing(){
    /*
        Descobre número de faces. (Somente para triângulos, por enquanto.)
    */
    set<pair<int,int>> faces; 
    for(auto it = this->elements.begin(); it != this->elements.end(); ++it){
        if(it->second.getElementType() == 2){
            
            vector<int>* nodes =  it->second.getNodes();
            
            /*Verificação das faces únicas*/
            int a = (*nodes)[0];
            int b = (*nodes)[1];
            int c = (*nodes)[2];
            
            //cria os pares de faces do elemento triangular
            //ordena sempre minmax para não ter problema de inserir (1,2) e (2,1) e contar como único no set, mas na verdade não ser único.
            pair<int, int> ab = minmax(a, b);
            pair<int, int> bc = minmax(b, c);
            pair<int, int> ca = minmax(c, a);

            faces.insert(ab);
            faces.insert(bc);
            faces.insert(ca);
            /*Fim verificação das faces únicas*/
            
            /*Cálculo do centroide do elemento*/
            Node n1 = this->nodes[a];
            Node n2 = this->nodes[b];
            Node n3 = this->nodes[c];

            double xc = (n1.getX() + n2.getX() + n3.getX()) / 3.0;
            double yc = (n1.getY() + n2.getY() + n3.getY()) / 3.0;
            double zc = (n1.getZ() + n2.getZ() + n3.getZ()) / 3.0;

            Node centroidNode = Node(-1, xc, yc, zc);
            this->centroids.emplace(it->first, centroidNode);
            /*Fim cálculo do centroide do elemento*/
        }
    }
    this->nfaces = faces.size();
    for(const auto& f : faces){
        this->faces.push_back(f);
        Node n1 = this->nodes[f.first];
        Node n2 = this->nodes[f.second];
        Node center = Node(-1, (n1.getX() + n2.getX())/2.0, (n1.getY() + n2.getY())/2.0, (n1.getZ() + n2.getZ())/2.0);
        this->faceMiddlePoints.push_back(center);
        double area = sqrt(pow(n1.getX() - n2.getX(), 2) + pow(n1.getY() - n2.getY(), 2) + pow(n1.getZ() - n2.getZ(), 2));
        this->faceAreas.push_back(area);
    }

    for(auto it = this->elements.begin(); it != this->elements.end(); ++it){
        if(it->second.getElementType() == 2){
            vector<int>* nodes =  it->second.getNodes();
            int a = (*nodes)[0];
            int b = (*nodes)[1];
            int c = (*nodes)[2];
            vector<pair<int,int>> edges;
            edges.push_back(minmax(a, b));
            edges.push_back(minmax(b, c));
            edges.push_back(minmax(c, a));
            for(int j = 0; j < edges.size(); j++){
                auto ite = find(faces.begin(), faces.end(), edges[j]);
                if(ite != faces.end()){
                    int index = distance(faces.begin(), ite);   
                    it->second.insertFace(index);
                }else{
                    cout << "ERROR: index not found!" << endl;
                }
            }
        }
    }
}

void Mesh::meshSummary(){
    cout << "Geometria do problema: " << this->geom_type << "D" << endl;
    cout << "Número de celulas (Contando triângulos, tetraedros e quadriláteros): " << this->ncells << endl;
    cout << "Número de faces (Somente funciona para triângulos): " << this->nfaces << endl;
    cout << "#-----------------------------------------------------------------------------------#" << endl;
    for(int i = 0; i < this->faces.size(); i++){
        cout << "Nós que compõem a face: " << this->faces[i].first << " " << this->faces[i].second << " | coordenadas do centro da face: " << this->faceMiddlePoints[i].getX() << " " << this->faceMiddlePoints[i].getY() << " " << this->faceMiddlePoints[i].getZ() << " | comprimento (2D por enquanto): " << this->faceAreas[i] << endl; 
    }
    
    cout << "#-----------------------------------------------------------------------------------#" << endl;
    cout << "Quantidade de nós: " << nnodes << endl;
    for(auto it = this->nodes.begin(); it != this->nodes.end(); ++it){
        cout << "Id do nó: " << it->first << "\tx: " << it->second.getX() << "\ty: " << it->second.getY() << "\tz: " << it->second.getZ() << endl; 
    }
    
    cout << "#-----------------------------------------------------------------------------------#" << endl;
    cout << "Quantidade total de elementos: " << totalElements << endl;
    for(auto it = this->elements.begin(); it != this->elements.end(); ++it){
        cout << "Id do elemento: " << it->first << "\ttipo do elemento: " << it->second.getElementType() << "\tid de seus nós: ";
        vector<int>* nodes = it->second.getNodes();
        for(int i = 0; i < nodes->size(); i++){
            cout << (*nodes)[i] << " ";
        }
        cout << endl;
    }    
    
    cout << "#-----------------------------------------------------------------------------------#" << endl;
    cout << "Quantidade total de grupos físicos: " << totalPhysicalEntities << endl;
    for(auto it = this->physicalGroups.begin(); it != this->physicalGroups.end(); ++it){
        cout << "Id do grupo: " << it->first << "\tdimensão: " << it->second.getDimension() << "\tnome: " << it->second.getName() << "\tid de seus elementos: ";
        vector<int>* elements = it->second.getElementIds();
        for(int i = 0; i < elements->size(); i++){
            cout << (*elements)[i] << " ";
        }
        cout << endl;
    }   

    cout << "#-----------------------------------------------------------------------------------#" << endl;
    cout << "Informação dos centroides:" << endl;
    for(auto it = this->centroids.begin(); it != this->centroids.end(); ++it){
        cout << "Id do elemento: " << it->first << "\t\txc: " << it->second.getX() << "\tyc: " << it->second.getY() << "\t\tzc: " << it->second.getZ() << endl; 
    }   
}