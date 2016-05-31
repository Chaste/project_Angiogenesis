import vtk
import vtk_tools.converter
from IPython.display import Image

class Scene():
    def __init__(self):
        self.renderer = vtk.vtkRenderer()
        self.renderer.SetBackground(1.0, 1.0, 1.0)
        
    def add_axes(self, polydata):
        axes = vtk.vtkCubeAxesActor2D()
        normals = vtk.vtkPolyDataNormals()
        normals.SetInput(polydata)
        axes.SetInput(normals.GetOutput())
        axes.SetCamera(self.renderer.GetActiveCamera())
        axes.SetLabelFormat("%6.4g")
        axes.SetFlyModeToOuterEdges()
        axes.SetFontFactor(2.0)
        tprop = vtk.vtkTextProperty()
        tprop.SetColor(1, 1, 1)
        axes.SetAxisTitleTextProperty(tprop)
        axes.SetAxisLabelTextProperty(tprop)
        self.renderer.AddViewProp(axes)
        
    def add_part(self, part, withAxes = True, dimension = 3):
        polydata = part.GetVtk(True)
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInput(polydata)
        actor = vtk.vtkActor()
        actor.GetProperty().SetColor(0,0,1) # (R,G,B)
        actor.SetMapper(mapper)
        self.renderer.AddActor(actor) 
        
    def add_points(self, points):
        converter = vtk_tools.converter.PointConverter()
        polydata = converter.convert(points)
        
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
        actor.GetProperty().SetColor(0,0,0) # (R,G,B)
        actor.SetMapper(mapper)
        
        self.renderer.AddActor(actor)
        
    def add_cells(self, cell_population):
        converter = vtk_tools.converter.CellConverter()
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
        actor.SetMapper(mapper)
        
        self.renderer.AddActor(actor)
                    
    def add_mesh(self, mesh, dimension = 3):
        converter = vtk_tools.converter.MeshConverter()
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
        
        self.renderer.AddActor(actor) 
        self.renderer.AddActor(actor2)  
        return actor, actor2
    
    def add(self, component):
        if "list" in component.__class__.__name__:
            component_type = type(component[0]).__name__
        else:
            component_type = type(component).__name__
        
        # Render the component, depending on its type
            if "Part" in component_type:
                self.add_part(component, False)
            
            if "Vertex" in component_type:
                self.add_points(component)
                
            if "SimpleCellPopulation" in component_type:
                self.sadd_cells(component)
            
            if "Mesh" in component_type:
                self.add_mesh(component)
    
    def get_renderer(self):
        return self.renderer
    
    def show(self, interactive = False, width=400, height=300):
        
        if interactive:
            self.renderer.ResetCamera()
            renWin = vtk.vtkRenderWindow()
            renWin.AddRenderer(self.renderer)
            renWin.SetSize(width, height)
            iren = vtk.vtkRenderWindowInteractor()
            iren.SetRenderWindow(renWin)
            iren.Initialize()
            iren.Start()
        else:
            renderWindow = vtk.vtkRenderWindow()
            renderWindow.SetOffScreenRendering(1)
            renderWindow.AddRenderer(self.renderer)
            renderWindow.SetSize(width, height)
            renderWindow.Render()
             
            windowToImageFilter = vtk.vtkWindowToImageFilter()
            windowToImageFilter.SetInput(renderWindow)
            windowToImageFilter.Update()
             
            writer = vtk.vtkPNGWriter()
            writer.SetWriteToMemory(1)
            writer.SetInputConnection(windowToImageFilter.GetOutputPort())
            writer.Write()
            data = str(buffer(writer.GetResult()))
            return Image(data)            