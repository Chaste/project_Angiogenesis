import vtk

class NetworkToPlanarSurface():
    
    """
    Take in a vessel network and return a VTK surface formed by the intersection with a cutting plane, assuming cylindrical vessels 
    
    The cutting plane can be described by extents (x, y), origin (x,y,z) and normal (x,y,z).
    
    This class uses FreeCAD to do the intersection with the plane
    """
    
    def __init__(self):
        
        self.Base = __import__("FreeCAD")
        self.Part = __import__("Part")
        self.network = None
        self.surface = None
        self.plane_extents = [1.e6,  1.e6]
        self.plane_origin = [-5.e-5, -5.e-5, 0.0]
        self.plane_normal = [0.0, 0.0, 1.0]
        
    def set_network(self, network):
        self.network = network
        
    def set_plane_extents(self, extents):
        self.plane_extents = extents
        
    def set_plane_normal(self, normal):
        self.plane_normal = normal
        
    def set_plane_origin(self, origin):
        self.plane_origin = origin
        
    def update(self):
        
        vessels = self.network.GetVessels()
        for eachVessel in vessels:
            
            start = eachVessel.GetStartNode().GetLocation()
            end = eachVessel.GetEndNode().GetLocation()
        
            cylinder = self.Part.makeCylinder(eachVessel.GetStartNode().GetRadius(), eachVessel.GetLength(), 
                                         self.Base.Vector(tuple(start)), self.Base.Vector(tuple(end - start)))
            cylinder = self.Part.makeShell(cylinder.Faces)
            
        clipping_plane = self.Part.makePlane(self.plane_extents[0], self.plane_extents[1], 
                                             self.Base.Vector(self.plane_origin[0], self.plane_origin[1], self.plane_origin[2]), 
                                             self.Base.Vector(self.plane_normal[0], self.plane_normal[1], self.plane_normal[2]))
        section = cylinder.section(clipping_plane)
        
        polyData = vtk.vtkPolyData() 
        points = vtk.vtkPoints()
        lines = vtk.vtkCellArray()
        for eachEdge in section.Edges:
            line = vtk.vtkLine()
            pointId1 = points.InsertNextPoint(eachEdge.Vertexes[0].Point)
            pointId2 = points.InsertNextPoint(eachEdge.Vertexes[1].Point)
            line.GetPointIds().InsertId(0, pointId1)
            line.GetPointIds().InsertId(1, pointId2)
            lines.InsertNextCell(line)
             
        polyData.SetPoints(points)
        polyData.SetLines(lines)   
        clean = vtk.vtkCleanPolyData()
        clean.SetInput(polyData)
        clean.Update()
        cylinder = clean.GetOutput()
        self.geometry = cylinder
            
        return self.geometry
    
    def get_output(self):
        
        return self.geometry
    
class NetworkTo3dCAD():
    
    """
    Take in a vessel network and return a 3d CAD representation, assuming cylindrical vessels
    This class uses FreeCAD to construct the CAD model
    """
    
    def __init__(self):
        
        self.Base = __import__("FreeCAD")
        self.Part = __import__("Part")
        self.network = None
        self.geometry = None
        
    def set_network(self, network):
        self.network = network
        
    def update(self):
        vessels = self.network.GetVessels()
        for eachVessel in vessels:
            
            start = eachVessel.GetStartNode().GetLocation()
            end = eachVessel.GetEndNode().GetLocation()
        
            cylinder = self.Part.makeCylinder(eachVessel.GetStartNode().GetRadius(), eachVessel.GetLength(), 
                                         self.Base.Vector(tuple(start)), self.Base.Vector(tuple(end - start)))
            cylinder = self.Part.makeShell(cylinder.Faces)
            
        self.geometry = cylinder
        return self.geometry
    
    def get_output(self):
        return self.geometry
    
class NetworkToVtkLines():
    
    def __init__(self):
        self.network = None
        self.surface = None
        
    def set_network(self, network):
        self.network = network
        
    def update(self):
        polyData = vtk.vtkPolyData() 
        points = vtk.vtkPoints()
        lines = vtk.vtkCellArray()
        info = vtk.vtkFloatArray()
        info.SetNumberOfComponents(1)
        info.SetName("Radius")
        
        for eachVessel in self.network.GetVessels():
            line = vtk.vtkLine()
            avRadius = 0
            for idx, eachSegment in enumerate(eachVessel.GetSegments()):
                pointId = points.InsertNextPoint(eachSegment.GetNodes().first.GetLocation())
                line.GetPointIds().InsertId(idx, pointId)
                avRadius = avRadius + eachSegment.GetRadius()
                if idx == len(eachVessel.GetSegments()) - 1:
                    pointId = points.InsertNextPoint(eachSegment.GetNodes().second.GetLocation())
                    line.GetPointIds().InsertId(idx + 1, pointId)
            avRadius = avRadius / len(eachVessel.GetSegments())
            lines.InsertNextCell(line)
            info.InsertNextTupleValue([avRadius])   
                
        polyData.SetPoints(points)
        polyData.SetLines(lines)
        polyData.GetCellData().SetScalars(info)
        self.surface = polyData
        return self.surface
    
    def get_output(self):
        
        return self.surface