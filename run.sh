#!/bin/bash

# caso dê erro.
set -e

# roda o solver
cd build
make
./TCC.exe

# volta para plotar o gráfico
cd ..
../ParaView-6.0.0-MPI-Linux-Python3.12-x86_64/bin/pvpython python/main.py