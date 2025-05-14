# TCC-UFJF
Repositório para o TCC.

Configuração das pastas:
<ul>
<li><b>app</b>: Arquivos principais vão aqui.</li>
<li><b>build</b>: Contém arquivos objetos, e é limpado usando <i>clean</i>.</li>
<li><b>docs</b>: Contém notas e outros arquivos para ajudar no desenvolvimento do compilador.</li>
<li><b>include</b>: Arquivos de cabeçalho.</li>
<li><b>inputs</b>: Arquivos fonte da linguagem para serem compilados e executados.</li>
<li><b>src</b>: Arquivo fonte que implementam os includes.</li>
</ul>

Guia básico de compilação:

```bash
# Na raiz do projeto. -S indica onde achar o arquivo fonte e -B indica onde fará o build.
$ cmake -S. -B ./build

# Verificando se o build deu certo, o output deverá ser algo como:
$ ls -l ./build
-rw-rw-r - 1 CMakeCache.txt
drwxrwxr-x 5 CMakeFiles
-rw-rw-r - 1 cmake_install.cmake
-rw-rw-r - 1 Makefile

# Agora entre no diretório build e dê make:
$ cd ./build
$ make
[ 50%] Building CXX object CMakeFiles/hello_world.dir/hello_world.cpp.o
[100%] Linking CXX executable hello_world
[100%] Built target hello_world

$ ./Compiladores.exe
Hello, World!
```