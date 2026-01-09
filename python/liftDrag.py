#!/usr/bin/env python3
import sys
import numpy as np
from collections import defaultdict

# ================== PARÂMETROS DO BENCHMARK ==================
CYL_CENTER = np.array([0.2, 0.2])
CYL_RADIUS = 0.05

RHO   = 1.0
U_INF = -0.2
RE    = 20.0

D  = 2.0 * CYL_RADIUS
NU = 1e-3      # viscosidade cinemática (Schäfer)
# =============================================================


# ------------------------------------------------------------
# Leitura simples de VTK ASCII (UNSTRUCTURED_GRID)
# ------------------------------------------------------------
def read_vtk(filename):
    with open(filename) as f:
        lines = f.readlines()

    points = []
    cells = []
    velocities = []
    pressures = []

    i = 0
    while i < len(lines):
        line = lines[i].strip()

        if line.startswith("POINTS"):
            n = int(line.split()[1])
            i += 1
            for _ in range(n):
                x, y, *_ = map(float, lines[i].split())
                points.append([x, y])
                i += 1

        elif line.startswith("CELLS"):
            n = int(line.split()[1])
            i += 1
            for _ in range(n):
                data = list(map(int, lines[i].split()))
                cells.append(data[1:])
                i += 1

        elif line.startswith("CELL_DATA"):
            i += 1
            while i < len(lines):
                l = lines[i].strip()

                if l.startswith("VECTORS"):
                    i += 1
                    for _ in range(len(cells)):
                        u, v, *_ = map(float, lines[i].split())
                        velocities.append([u, v])
                        i += 1

                elif l.startswith("SCALARS"):
                    i += 2
                    for _ in range(len(cells)):
                        pressures.append(float(lines[i]))
                        i += 1
                else:
                    i += 1
        else:
            i += 1

    return (
        np.array(points),
        cells,
        np.array(velocities),
        np.array(pressures),
    )


# ------------------------------------------------------------
# Geometria auxiliar
# ------------------------------------------------------------
def centroid(cell, points):
    return np.mean(points[cell], axis=0)


def build_edge_map(cells):
    edge_map = defaultdict(list)
    for ci, cell in enumerate(cells):
        n = len(cell)
        for k in range(n):
            a, b = cell[k], cell[(k + 1) % n]
            edge = tuple(sorted((a, b)))
            edge_map[edge].append((ci, a, b))
    return edge_map


def boundary_edges(edge_map):
    return {e: v[0] for e, v in edge_map.items() if len(v) == 1}


# ------------------------------------------------------------
# Gradientes por least squares (centro–centro)
# ------------------------------------------------------------
def compute_gradients(centers, values, neighbors):
    grads = np.zeros((len(centers), 2))
    for i, nbrs in enumerate(neighbors):
        if len(nbrs) < 2:
            continue
        A = []
        b = []
        for j in nbrs:
            dx = centers[j] - centers[i]
            A.append(dx)
            b.append(values[j] - values[i])
        A = np.array(A)
        b = np.array(b)
        grads[i], *_ = np.linalg.lstsq(A, b, rcond=None)
    return grads


# ------------------------------------------------------------
# Integração das forças (Schäfer & Turek)
# ------------------------------------------------------------
def compute_forces(points, cells, velocities, pressures):
    centers = np.array([centroid(c, points) for c in cells])

    edge_map = build_edge_map(cells)
    bnd_edges = boundary_edges(edge_map)

    neighbors = [[] for _ in cells]
    for owners in edge_map.values():
        if len(owners) == 2:
            c1, c2 = owners[0][0], owners[1][0]
            neighbors[c1].append(c2)
            neighbors[c2].append(c1)

    grad_u = compute_gradients(centers, velocities[:, 0], neighbors)
    grad_v = compute_gradients(centers, velocities[:, 1], neighbors)

    Fp = np.zeros(2)
    Fv = np.zeros(2)

    for (a, b), (ci, _, _) in bnd_edges.items():
        pa, pb = points[a], points[b]
        mid = 0.5 * (pa + pb)

        # verifica se pertence ao cilindro
        r = np.linalg.norm(mid - CYL_CENTER)
        if abs(r - CYL_RADIUS) > 1e-3:
            continue

        edge = pb - pa
        ds = np.linalg.norm(edge)
        if ds == 0.0:
            continue

        # normal externa do cilindro
        n = mid - CYL_CENTER
        n /= np.linalg.norm(n)

        # ---------------- PRESSÃO ----------------
        p_face = pressures[ci]
        Fp += -p_face * n * ds

        # ---------------- VISCOUS (2νD(u)) -------
        grad = np.array([
            [grad_u[ci, 0], grad_u[ci, 1]],
            [grad_v[ci, 0], grad_v[ci, 1]]
        ])

        Dmat = 0.5 * (grad + grad.T)
        tau = 2.0 * NU * Dmat

        Fv += tau @ n * ds

    return Fp + Fv, Fp, Fv


# ------------------------------------------------------------
# MAIN
# ------------------------------------------------------------
def main(filename):
    points, cells, velocities, pressures = read_vtk(filename)

    F, Fp, Fv = compute_forces(points, cells, velocities, pressures)

    Cd = 2.0 * F[0] / (RHO * U_INF**2 * D)
    Cl = 2.0 * F[1] / (RHO * U_INF**2 * D)

    print("=== BENCHMARK Schäfer & Turek – 2D-1 (Re = 20) ===")
    print(f"Força total   : Fx = {F[0]:.6e}, Fy = {F[1]:.6e}")
    print(f"Força pressão : Fx = {Fp[0]:.6e}, Fy = {Fp[1]:.6e}")
    print(f"Força viscosa : Fx = {Fv[0]:.6e}, Fy = {Fv[1]:.6e}")
    print(f"Cd = {Cd:.6f}")
    print(f"Cl = {Cl:.6f}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Uso: python compute_cd_cl_schafer.py arquivo.vtk")
        sys.exit(1)
    main(sys.argv[1])
