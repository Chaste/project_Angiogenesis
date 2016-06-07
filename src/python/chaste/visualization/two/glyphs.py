"""
2D Glyphs which can be attached to matploblib plots and interacted with
"""

import matplotlib.patches

import chaste.geometry.converters.other

class VertexGlyph():
    
    def __init__(self, vertex):
        self.vertex = vertex
        self.patch = matplotlib.patches.Circle(self.vertex.get_location(), 
                                               radius = 0.2, 
                                               color = 'blue', 
                                               alpha = 0.75, 
                                               picker=True)
        
    def attach_to_axes(self, axes):
        axes.add_patch(self.patch)
        
class PolygonGlyph():
    
    def __init__(self, polygon):
        self.polygon = polygon
        x = []
        y = []
        vertices = polygon.get_vertices()
        for eachVertex in vertices:
            x.append(eachVertex.get_location()[0])
            y.append(eachVertex.get_location()[1])
        self.patch = matplotlib.patches.Polygon([x, y], color = 'red',  alpha = 0.75, picker = True)
        
    def attach_to_axes(self, axes):
        axes.add_patch(self.patch)
        
class PartGlyph():
    
    def __init__(self, part):
        polygons = part.get_polygons()
        self.patches = []
        for eachPolygon in polygons:
            x = []
            y = []
            vertices = eachPolygon.get_vertices()
            for eachVertex in vertices:
                x.append(eachVertex.get_location()[0])
                y.append(eachVertex.get_location()[1])
            self.patches.append(matplotlib.patches.Polygon([x, y], color = 'red',  alpha = 0.75, picker = True))
        
    def attach_to_axes(self, axes):
        for eachPatch in self.patches:
            axes.add_patch(eachPatch)
            
class MeshGlyph():
    
    def __init__(self, mesh, chaste_format = True, color = "blue"):
        
        if chaste_format:
            node_locations = mesh.get_node_locations()
            connectivity = mesh.get_connectivity()
        else:
            node_locations = mesh[0]
            connectivity = mesh[1]
        
        lines = []
        for eachElement in connectivity:
            for idx in range(len(eachElement)-1):
                lines.append((node_locations[eachElement[idx]], node_locations[eachElement[idx+1]]))
            lines.append((node_locations[eachElement[-1]], node_locations[eachElement[0]]))
        self.pacthes = matplotlib.collections.LineCollection(lines, color=color) 
        
    def attach_to_axes(self, axes):
        axes.add_collection(self.pacthes) 
        
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
        
class VtkLinesGlyph():
    
    def __init__(self, surface, color="green"):
        
        converter = chaste.geometry.other.VtkToTri()
        boundary = converter.generate(surface)
        
        node_locations = boundary[0]
        connectivity = boundary[1]
        
        lines = []
        for eachElement in connectivity:
            for idx in range(len(eachElement)-1):
                lines.append((node_locations[eachElement[idx]], node_locations[eachElement[idx+1]]))
            lines.append((node_locations[eachElement[-1]], node_locations[eachElement[0]]))
        self.pacthes = matplotlib.collections.LineCollection(lines, color=color) 
        
    def attach_to_axes(self, axes):
        axes.add_collection(self.pacthes) 
        
class CellGlyph():
    def __init__(self, cell):
        self.cell = cell
        self.position = cell.get_location()
        self.radius = 5.0
        self.inner_circle = matplotlib.patches.Circle(self.position, radius = self.radius/4.0, color="yellow", alpha = 0.8)
        self.circle = matplotlib.patches.Circle(self.position, radius = self.radius, color="green", alpha = 0.5, picker=True)
        
    def attach_to_axes(self, axes):
        axes.add_patch(self.circle)
        axes.add_patch(self.inner_circle)
        
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