from __future__ import print_function
import sys

if(len(sys.argv) != 3):
  print('Usage: ', sys.argv[0], 'input.obj output.vtk')
  sys.exit()

import vtk
reader = vtk.vtkOBJReader()
reader.SetFileName(sys.argv[1])
reader.Update()
obj = reader.GetOutput()

writer = vtk.vtkPolyDataWriter()
writer.SetFileName(sys.argv[2])
#writer.SetInput(obj)
#if vtk.VTK_MAJOR_VERSION <= 100:
    
#else:
writer.SetInputData(obj)
writer.Write()
