import matplotlib.patches
import vtk

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
        
class CellGlyph3d():
    
    def __init__(self, cell_population):
        converter = casie.rwc.vtk_converter.CellConverter()
        polydata = converter.convert(cell_population)
        
        # Each point is a sphere
        sphere = vtk.vtkSphereSource()
        sphere.SetRadius(1.0)
        sphere.SetPhiResolution(16)
        sphere.SetThetaResolution(16)
         
        # make the glyphs
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
        actor.GetProperty().SetColor(0,1,0) # (R,G,B)
        self.actor = actor