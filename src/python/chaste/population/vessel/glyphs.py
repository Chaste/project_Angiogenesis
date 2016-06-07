import matplotlib.patches
import vtk

class VesselNetworkGlyph():
    
    def __init__(self, network):
        
        node_locations = []
        nodes = network.nodes
        for idx, eachNode in enumerate(nodes):
            eachNode.id = idx
            node_locations.append((eachNode.GetLocationVector()[0], eachNode.GetLocationVector()[1]))
            
        connectivity = []
        for eachVessel in network.vessels:
            connectivity.append((eachVessel.start_node.id, eachVessel.end_node.id))
        
        lines = []
        for eachElement in connectivity:
            for idx in range(len(eachElement)-1):
                lines.append((node_locations[eachElement[idx]], node_locations[eachElement[idx+1]]))
            lines.append((node_locations[eachElement[-1]], node_locations[eachElement[0]]))
        self.pacthes = matplotlib.collections.LineCollection(lines, color='red') 
        
    def attach_to_axes(self, axes):
        axes.add_collection(self.pacthes) 
        
class NodeGlyph():
    def __init__(self, node):
        self.node = node
        self.position = node.get_location()
        self.color = 'red'
        self.radius = 2.0
        self.circle = matplotlib.patches.Circle(self.position, radius = self.radius, color=self.color, alpha = 0.75, picker=True)
        
    def attach_to_axes(self, axes):
        axes.add_patch(self.circle)
        
class SegmentGlyph():
    def __init__(self, segment):
        self.segment = segment
        self.start_position = segment.startNode.location
        self.end_position = segment.endNode.location
        self.color = 'red'
        self.radius = segment.radius
        self.line = matplotlib.lines.Line2D([self.start_position[0], self.end_position[0]], 
                                     [self.start_position[1], self.end_position[1]], 
                                     picker=True, 
                                     color = 'red', alpha = 0.5, lw = self.radius)
        
    def attach_to_axes(self, axes):
        axes.add_line(self.line)
        
class VesselGlyph3d():
    
    def __init__(self, vessel_network):
        pointColors = vtk.vtkUnsignedCharArray()
        pointColors.SetNumberOfComponents(3)
        pointColors.SetName("PointColors") 
         
        vessel_vtk = vessel_network.GetVtk()
        for _ in range(vessel_vtk.GetNumberOfPoints()):
            pointColors.InsertNextTupleValue([1, 1, 1])  
         
        sphere = vtk.vtkSphereSource()
        sphere.SetRadius(0.5)
        sphere.SetPhiResolution(16)
        sphere.SetThetaResolution(16)
         
        glyph = vtk.vtkGlyph3D()
        glyph.SetInput(vessel_vtk)
        glyph.SetSource(sphere.GetOutput())
        glyph.ClampingOff()
        glyph.SetScaleModeToScaleByScalar()
        glyph.SetScaleFactor(1.0) 
         
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInput(glyph.GetOutput())
        mapper.ScalarVisibilityOn()
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)        
         
        colors = vtk.vtkUnsignedCharArray()
        colors.SetNumberOfComponents(3)
        colors.SetName("Colors") 
        for _ in range(vessel_vtk.GetNumberOfCells()):
            colors.InsertNextTupleValue([1, 1, 1])      
        vessel_vtk.GetCellData().SetScalars(colors) 
        
        mapper2 = vtk.vtkPolyDataMapper()
        mapper2.SetInput(vessel_vtk)
        mapper2.ScalarVisibilityOn()
        actor2 = vtk.vtkActor()
        actor2.SetMapper(mapper2)
        self.actors = [actor, actor2]