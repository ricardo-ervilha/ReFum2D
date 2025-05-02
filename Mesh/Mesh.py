import gmsh
import sys

# initialize gmsh
gmsh.initialize()

# square points
lc = 1e-1
p1 = gmsh.model.geo.add_point(0, 0, 0, lc)
p2 = gmsh.model.geo.add_point(0, 1, 0, lc)
p3 = gmsh.model.geo.add_point(1, 1, 0, lc)
p4 = gmsh.model.geo.add_point(1, 0, 0, lc)

# edges
l1 = gmsh.model.geo.add_line(p1, p2)
l2 = gmsh.model.geo.add_line(p2, p3)
l3 = gmsh.model.geo.add_line(p3, p4)
l4 = gmsh.model.geo.add_line(p4, p1)

# create face
loop = gmsh.model.geo.add_curve_loop([l1, l2, l3, l4])
gmsh.model.geo.add_plane_surface([loop])

gmsh.model.geo.synchronize()
gmsh.model.mesh.generate(2)

gmsh.write("GFG.msh")

# finalize gmsh
gmsh.finalize()
