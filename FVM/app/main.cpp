#include <gmsh.h>
#include <iostream>
#include <vector>

using namespace std;

int main(int argc, char **argv) {
    string filename = "../inputs/mesh.msh";
    
    /*Inicia Leitura do Arquivo*/
    gmsh::initialize();
    try{
        gmsh::open(filename); //tenta abrir o arquivo
    }catch (const std::exception &e) {
        std::cerr << "Erro ao abrir o arquivo: " << e.what() << std::endl;
        gmsh::finalize();
        exit(1);
    }

    std::vector<std::size_t> nodeTags;
    std::vector<double> nodeCoords;
    std::vector<double> nodeParametricCoords;

    gmsh::model::mesh::getNodes(nodeTags, nodeCoords, nodeParametricCoords);

    std::cout << "Quantidade de nós: " << nodeTags.size() << "\n";

    for (size_t i = 0; i < nodeTags.size(); ++i) {
        std::cout << "Nó " << nodeTags[i] << ": ("
                    << nodeCoords[3 * i] << ", "
                    << nodeCoords[3 * i + 1] << ", "
                    << nodeCoords[3 * i + 2] << ")\n";
    }

    std::vector<int> elementTypes;
    std::vector<std::vector<std::size_t>> elementTags, nodeTagsPerElem;

    gmsh::model::mesh::getElements(elementTypes, elementTags, nodeTagsPerElem);

    std::cout << "Número de tipos de elementos: " << elementTypes.size() << "\n";

    for (size_t i = 0; i < elementTypes.size(); ++i) {
        std::cout << "Tipo de elemento " << elementTypes[i]
                    << " com " << elementTags[i].size()
                    << " elementos\n";
    }

    gmsh::finalize();

    return 0;
}
