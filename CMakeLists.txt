cmake_minimum_required(VERSION 3.10)
project(TCC)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})

set(CMAKE_CXX_FLAGS_DEBUG "-pg -O3 -g")

include_directories(include)

# Fontes e cabeçalhos
set(SOURCE_FILES
    src/Mesh.cpp
    # src/FVMSolver.cpp
    # src/Diffusion.cpp
    # src/Source.cpp
    # src/GradientReconstruction.cpp
    # src/Convection.cpp
    src/NSSolver.cpp
)
set(HEADER_FILES
    include/Element.h
    include/PhysicalEntity.h
    include/Node.h
    include/Edge.h
    include/Mesh.h
    include/Cell.h  
    # include/FVMsolver.h
    include/utils.h
    include/BoundaryCondition.h
    # include/Diffusion.h
    # include/Source.h
    # include/GradientReconstruction.h
    # include/Convection.h
    include/NSSolver.h
)

add_executable(${PROJECT_NAME}.exe app/main.cpp ${SOURCE_FILES})
target_link_libraries(${PROJECT_NAME}.exe ${ARMADILLO_LIBRARIES})
target_precompile_headers(${PROJECT_NAME}.exe PRIVATE include/pch.h)