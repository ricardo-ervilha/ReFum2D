# TCC - UFJF

**Aluno:** Ricardo Ervilha Silva  
**Orientador:** Prof. José Jerônimo Camata

_Repositório do Trabalho de Conclusão de Curso (2025.3) da Universidade Federal de Juiz de Fora - Ciência da Computação_

![Status: Em desenvolvimento](https://img.shields.io/badge/status-Em%20Desenvolvimento-yellow)
![GMSH v2](https://img.shields.io/badge/gmsh-v2%20ASCII-blue)

---

## 🎯 Objetivo

Este projeto implementa um solucionador numérico baseado no **Método dos Volumes Finitos (FVM)**, focado na leitura e simulação de malhas `.msh` não estruturadas geradas pelo Gmsh (**versão 2 ASCII**).

As geometrias devem ser criadas com **pontos, linhas e condições de contorno na ordem anti-horária**, o que facilita o processamento e interpretação dos dados da malha.

---

## 📁 Estrutura de Diretórios
```
📁 Estrutura do Projeto

├── app/        # Contém o ponto de entrada principal do programa (main.cpp)
├── build/      # Diretório gerado automaticamente com os arquivos de build pelo CMake
├── docs/       # Documentação geral, imagens ilustrativas e arquivos auxiliares
├── include/    # Arquivos de cabeçalho (.h) com definições de classes e interfaces
├── inputs/     # Malhas de entrada no formato .msh (versão 2 ASCII, geradas com Gmsh)
├── outputs/    # Arquivos de saída (.vtk) contendo os resultados da simulação para visualização no ParaView
├── src/        # Implementações das classes e funções declaradas nos headers (arquivos .cpp)
```

## ⚙️ Como utilizar o projeto

### Pré-requisitos

- G++ (ou outro compilador C++ compatível)
- CMake (>= 3.10)
- [Gmsh](https://gmsh.info/) para gerar as malhas `.msh`
- [ParaView](https://www.paraview.org/) (opcional, para visualização dos resultados `.vtk`)

### Passos

```bash
# Na raiz do projeto
mkdir build
cd build
cmake ..
make

# Após compilado com sucesso
./TCC.exe

../../ParaView-6.0.0-MPI-Linux-Python3.12-x86_64/bin/pvpython main.py