import subprocess
import os

build_dir = os.path.join("..", "build")
exe_path = os.path.join(build_dir, "TCC.exe")

inputs_dir = os.path.join("..", "inputs", "regular")
outputs_dir = os.path.join("..", "outputs")

meshes = [
    "quad_10x10.msh",
    "quad_20x20.msh",
    "quad_30x30.msh",
    "quad_40x40.msh",
    "quad_50x50.msh",
    "quad_100x100.msh",
    "quad_200x200.msh",
]

os.makedirs(outputs_dir, exist_ok=True)

for mesh in meshes:
    # Extrai sufixo (ex: "10x10" de "quad_10x10.msh")
    suffix = mesh.replace("quad_", "").replace(".msh", "")

    # Caminhos de entrada e saída
    msh_path = os.path.join(inputs_dir, mesh)
    vtk_path = os.path.join(outputs_dir, f"result_{suffix}.vtk")

    # argv[1] = msh_path, argv[2] = vtk_path, argv[3] = suffix
    cmd = [exe_path, msh_path, vtk_path, suffix]

    print(f"Rodando: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)
