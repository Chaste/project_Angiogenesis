from FreeCAD import Base
import FreeCAD as App
import Part
import Mesh

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
        w2 = Pa
        
        sides = w1.extrude(Base.Vector(0,0,depth))
        front = Part.Face(w1)
        back = Part.Face(w2)
        
        domain = Part.Compound([sides, front, back])
        domain.translate(Base.Vector(60.0, 300.0, 45.0))
        
        compound = domain.fuse(shape)
        #compound.exportBrep(working_directory + "/merge.brep")
        compound.exportStl(working_directory + "/merge.stl")
        