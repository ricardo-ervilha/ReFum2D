#!/usr/bin/env pvpython

from paraview.simple import *

# janela onde será renderizado
renderView = GetActiveViewOrCreate('RenderView')

# carrega o vtk
vtkFile = LegacyVTKReader(FileNames=['solution.vtk'])

# filtro
cellToPoint = CellDatatoPointData(Input=vtkFile)

display = Show(cellToPoint, renderView)

# Renderiza
Interact()
Render()