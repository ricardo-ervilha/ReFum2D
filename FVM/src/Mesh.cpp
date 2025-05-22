#include "Mesh.h"
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
                }
            }
        }
        file.close();
    }
    else{
        cout << "ERROR: could not open the file." << endl;
        exit(1);
    }
}

void Mesh::meshSummary(){
    cout << "#-----------------------------------------------------------------------------------#" << endl;
    cout << "Quantidade de nós: " << nnodes << endl;
    for(auto it = this->nodes.begin(); it != this->nodes.end(); ++it){
        cout << "Id do nó: " << it->first << "\tx: " << it->second.getX() << "\ty: " << it->second.getY() << "\tz: " << it->second.getZ() << endl; 
    }
    
    cout << "#-----------------------------------------------------------------------------------#" << endl;
    cout << "Quantidade total de elemntos: " << totalElements << endl;
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
}