import matplotlib.patches
import vtk

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
            
class VertexGlyph3d():
    
    def __init__(self, vertex):
        self.vertex = vertex
        
        converter = casie.rwc.vtk_converter.PointConverter()
        polydata = converter.convert([vertex])
        sphere = vtk.vtkSphereSource()
        sphere.SetRadius(1.0)
        sphere.SetPhiResolution(16)
        sphere.SetThetaResolution(16)
        glyph = vtk.vtkGlyph3D()
        glyph.SetInput(polydata)
        glyph.SetSource(sphere.GetOutput())
        glyph.ClampingOff()
        glyph.SetScaleModeToScaleByScalar()
        glyph.SetScaleFactor(1.0) 
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInput(glyph.GetOutput())
        mapper.ScalarVisibilityOn()
        actor = vtk.vtkActor()
        actor.GetProperty().SetColor(0,0,0)
        actor.SetMapper(mapper)
        
        self.actor = actor
        
class VerticesGlyph3d():
    
    def __init__(self, vertices):
        self.vertices = vertices
        
        converter = casie.rwc.vtk_converter.PointConverter()
        polydata = converter.convert(vertices)
        sphere = vtk.vtkSphereSource()
        sphere.SetRadius(1.0)
        sphere.SetPhiResolution(16)
        sphere.SetThetaResolution(16)
        glyph = vtk.vtkGlyph3D()
        glyph.SetInput(polydata)
        glyph.SetSource(sphere.GetOutput())
        glyph.ClampingOff()
        glyph.SetScaleModeToScaleByScalar()
        glyph.SetScaleFactor(1.0) 
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInput(glyph.GetOutput())
        mapper.ScalarVisibilityOn()
        actor = vtk.vtkActor()
        actor.GetProperty().SetColor(0,0,0)
        actor.SetMapper(mapper)
        self.actor = actor
        
class PartGlyph3d():
    
    def __init__(self, part):
        polydata = part.GetVtk(True)
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInput(polydata)
        actor = vtk.vtkActor()
        actor.GetProperty().SetColor(0,0,1)
        actor.SetMapper(mapper)
        self.actor = actor