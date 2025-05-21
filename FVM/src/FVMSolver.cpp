#include "FVMSolver.h"
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>

FVMSolver::FVMSolver(){

}

FVMSolver::~FVMSolver(){

}

void FVMSolver::readMesh(string filepath){
    ifstream mesh(filepath);

    string line;

    if(mesh.is_open()){
        while(getline(mesh, line)){
            line.erase(remove(line.begin(), line.end(), '\r'), line.end()); /*remove caracteres problemáticos*/
            if (line == "$MeshFormat") {
                getline(mesh, line);
                istringstream iss(line);
                double version;
                iss >> version;
                if(version == 2.2)
                    continue;
                else
                    cout << "ERROR: Unexpected version of .msh file." << endl;
            } else if(line == "$PhysicalNames"){
                getline(mesh,line);
                istringstream iss(line);
                int totalPhysicalEntities;
                iss >> totalPhysicalEntities;
            }
        }
        mesh.close();
    }
    else{
        cout << "ERROR: could not open the file." << endl;
        exit(1);
    }
}