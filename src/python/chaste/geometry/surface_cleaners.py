import vtk
from vmtk import pypes

class SurfaceCleaner3d():
    
    def __init__(self, surface_file, work_dir, clipping_planes = None):
        
        self.surface_file = surface_file
        self.clipping_planes = clipping_planes
        self.work_dir = work_dir
        
    def generate(self, decimate_factor = 0.5):
        
        # read surface
        reader = vtk.vtkXMLPolyDataReader()
        reader.SetFileName(self.surface_file)
        reader.Update()
        
        # clip surface
        previous = reader
        if self.clipping_planes is not None:
            for eachPlane in self.clipping_planes:
                plane = vtk.vtkPlane()
                plane.SetOrigin(eachPlane[0])
                plane.SetNormal(eachPlane[1])
                
                clipper = vtk.vtkClipPolyData()
                clipper.SetInputConnection(previous.GetOutputPort())
                clipper.SetClipFunction(plane)
                clipper.SetValue(0.0)
                
                previous= clipper
                 
        # decimate
        deci = vtk.vtkDecimatePro()
        deci.SetInputConnection(previous.GetOutputPort())
        deci.SetTargetReduction(decimate_factor)
        
        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetFileName(self.work_dir + "/clipped.vtp")
        writer.SetInputConnection(deci.GetOutputPort())
        writer.Update() 
        
        # mesh
        input_location = self.work_dir + "clipped.vtp"
        outoput_location = self.work_dir + "remesh.vtu"
        myArguments = "vmtkmeshgenerator -ifile " + input_location + " -ofile " + outoput_location + " -edgelength 5"
        myPype = pypes.PypeRun(myArguments) # expensive, uncomment to do remeshing
        
        # get mesh surface
        input_location = self.work_dir + "remesh.vtu"
        outoput_location = self.work_dir + "remesh_surface.vtp"
        myArguments = "vmtkmeshtosurface -ifile " + input_location + " -ofile " + outoput_location + " -cleanoutput 1"
        myPype = pypes.PypeRun(myArguments)    
        
        # convert to stl
        reader = vtk.vtkXMLPolyDataReader()
        reader.SetFileName(self.work_dir + "remesh_surface.vtp")
        reader.Update()
        
        # Open cap version for centrelines
        # clip surface
        previous = reader
        if self.clipping_planes is not None:
            for eachPlane in self.clipping_planes:
                plane = vtk.vtkPlane()
                plane.SetOrigin(eachPlane[0])
                plane.SetNormal(eachPlane[1])
                
                clipper = vtk.vtkClipPolyData()
                clipper.SetInputConnection(previous.GetOutputPort())
                clipper.SetClipFunction(plane)
                clipper.SetValue(0.0)
                
                previous= clipper
        
        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetFileName(self.work_dir + "cleaned.vtp")
        writer.SetInputConnection(previous.GetOutputPort())
        writer.Update()     
        
        writer = vtk.vtkSTLWriter()
        writer.SetFileName(self.work_dir + "cleaned.stl")
        writer.SetInputConnection(reader.GetOutputPort())
        writer.SetFileTypeToASCII()
        writer.Update()   