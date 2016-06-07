import matplotlib
import chaste.mesh.converters

class VtkLinesGlyph():
    
    def __init__(self, surface, color="green"):
        converter = chaste.mesh.converters.VtkToTriMesh()
        converter.input = surface
        boundary = converter.update()
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