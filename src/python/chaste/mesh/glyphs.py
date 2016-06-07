import vtk
import matplotlib

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

class MeshGlyph3d():
    
    def __init__(self, mesh, dimension = 3):
        converter = casie.rwc.vtk_converter.MeshConverter()
        mesh_vtk = converter.convert(mesh, dimension)
        
        geo_filter = vtk.vtkGeometryFilter()
        geo_filter.SetInput(mesh_vtk)
        
        extract_edges = vtk.vtkExtractEdges()
        extract_edges.SetInput(mesh_vtk)
        
        tube_filter = vtk.vtkTubeFilter()
        tube_filter.SetInput(extract_edges.GetOutput())
        tube_filter.SetRadius(0.02)
        
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInput(geo_filter.GetOutput())
        #mapper.ScalarVisibilityOn()
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().SetColor(1,0,0) # (R,G,B)
        
        mapper2 = vtk.vtkPolyDataMapper()
        mapper2.SetInput(tube_filter.GetOutput())
        #mapper2.ScalarVisibilityOn()
        actor2 = vtk.vtkActor()
        actor2.SetMapper(mapper2)
        actor.GetProperty().SetColor(1,1,1) # (R,G,B
        
        self.actors = [actor, actor2]