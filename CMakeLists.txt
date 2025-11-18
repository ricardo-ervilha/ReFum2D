cmake_minimum_required(VERSION 3.10)
project(TCC)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})

find_package(yaml-cpp REQUIRED)

set(CMAKE_CXX_FLAGS_DEBUG "-O3")

include_directories(include)

# Adiciona Eigen (header-only)
set(EIGEN3_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/eigen-5.0.0")
include_directories(${EIGEN3_INCLUDE_DIR})

# Fontes e cabeçalhos
set(SOURCE_FILES
    src/Mesh.cpp
    src/Orchestrator.cpp
    src/ReFumSolver.cpp
)
set(HEADER_FILES
    include/Element.h
    include/PhysicalEntity.h
    include/Node.h
    include/Edge.h
    include/Mesh.h
    include/Cell.h  
    include/Orchestrator.h
    include/utils.h
    include/BoundaryCondition.h
    include/ReFumSolver.h
)

add_executable(${PROJECT_NAME}.exe app/main.cpp ${SOURCE_FILES})
target_link_libraries(${PROJECT_NAME}.exe PRIVATE ${ARMADILLO_LIBRARIES} yaml-cpp)
target_precompile_headers(${PROJECT_NAME}.exe PRIVATE include/pch.h)

# Alternativa moderna (opcional)
# target_include_directories(${PROJECT_NAME}.exe PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/eigen-5.0.0")
