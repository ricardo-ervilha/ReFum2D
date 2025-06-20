#!/bin/bash

# Caminho do arquivo VTK gerado no WSL
ORIGEM="../outputs/result.vtk"

# Caminho de destino no Windows montado no WSL
DESTINO="/mnt/c/Users/ricar/OneDrive/Documentos/results"

# Cria o diretório se não existir
mkdir -p "$DESTINO"

# Copia o arquivo
cp "$ORIGEM" "$DESTINO"

# Mensagem de sucesso
echo "✅ Arquivo copiado para: $DESTINO/result.vtk"