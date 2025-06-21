cmake_minimum_required(VERSION 3.10)
project(TCC)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# diretório do eigen
include_directories(eigen-3.4.0)

include_directories(include)
include_directories(${GMSH_SDK_DIR}/include)

link_directories(${GMSH_SDK_DIR}/lib)

set(SOURCE_FILES
    src/Mesh.cpp
    src/FVMSolver.cpp
)

set(HEADER_FILES
    include/Element.h
    include/PhysicalGroup.h
    include/Node.h
    include/Mesh.h
    include/FVMSolver.h
    include/BoundaryCondition.h
)

add_executable(${PROJECT_NAME}.exe app/main.cpp ${SOURCE_FILES})