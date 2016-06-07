import vtk
from FreeCAD import Base
import FreeCAD as App
import Part
import Mesh
from vmtk import pypes

class Cleaner3d():
    
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
        
class MergeStl():
    
    def __init__(self):
        
        pass
    
    def generate(self, surface_file, file_name, working_directory, geometry_tolerance = 0.05):
        
        Mesh.open(surface_file)
    
        # Convert to solid
        shape = Part.Shape()
        doc=App.activeDocument() 
        shape.makeShapeFromMesh(doc.getObject(file_name).Mesh.Topology, geometry_tolerance) 
        #solid = Part.makeSolid(shape)
        
        # Add bounding domain
        width = 480.0
        height = 180.0
        depth = 160.0
        
        e1 = Part.makeLine((0,0,0), (width,0,0))
        e2 = Part.makeLine((width,0,0), (width,height,0))
        e3 = Part.makeLine((width,height,0), (0,height,0))
        e4 = Part.makeLine((0,height,0), (0,0,0))
        e5 = Part.makeLine((0,0,depth), (width,0,depth))
        e6 = Part.makeLine((width,0,depth), (width,height,depth))
        e7 = Part.makeLine((width,height,depth), (0,height,depth))
        e8 = Part.makeLine((0,height,depth), (0,0,depth))
        w1 = Part.Wire([e1,e2,e3,e4]) 
        w2 = Part.Wire([e5,e6,e7,e8]) 
        
        sides = w1.extrude(Base.Vector(0,0,depth))
        front = Part.Face(w1)
        back = Part.Face(w2)
        
        domain = Part.Compound([sides, front, back])
        domain.translate(Base.Vector(60.0, 300.0, 45.0))
        
        compound = domain.fuse(shape)
        #compound.exportBrep(working_directory + "/merge.brep")
        compound.exportStl(working_directory + "/merge.stl")
        