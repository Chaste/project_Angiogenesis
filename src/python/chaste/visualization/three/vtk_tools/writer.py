import vtk

from casie.vtk_tools import converter
 
class MeshWriter():
     
    def __init__(self):
        pass
     
    def write(self, mesh, fileName, dimension = 2):
        
        converter = converter.MeshConverter()
        data = converter.convert(mesh, dimension)
        
        # Write the grid to file   
        writer = vtk.vtkXMLUnstructuredGridWriter() 
        writer.SetFileName(fileName)
        if vtk.VTK_MAJOR_VERSION <= 5:
            writer.SetInput(data)
        else:
            writer.SetInputData(data)
        writer.Write()
        
class VesselWriter():
    
    def __init__(self):
        pass
    
    def write(self, vesselNetwork, fileName):
        converter = converter.VesselConverter()
        data = converter.convert(vesselNetwork)
        
        # Write the file
        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetFileName(fileName)
        if vtk.VTK_MAJOR_VERSION <= 5:
            writer.SetInput(data)
        else:
            writer.SetInputData(data)
        writer.Write() 
        
class CellWriter():
    
    def __init__(self):
        pass
    
    def write(self, cellPopulation, fileName):
        converter = converter.CellConverter()
        data = converter.convert(cellPopulation)
        
        # Write the file
        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetFileName(fileName)
        if vtk.VTK_MAJOR_VERSION <= 5:
            writer.SetInput(data)
        else:
            writer.SetInputData(data)
        writer.Write() 
        
class PartWriter():
    
    def __init__(self):
        pass
    
    def write(self, part, fileName):
        data = part.get_vtk()
        
        # Write the file
        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetFileName(fileName)
        if vtk.VTK_MAJOR_VERSION <= 5:
            writer.SetInput(data)
        else:
            writer.SetInputData(data)
        writer.Write()   