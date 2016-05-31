"""
3D Glyphs which can be attached to vtk renderers and interacted with
"""

import vtk
import casie.rwc.vtk_converter

class VertexGlyph():
    
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
        
class VerticesGlyph():
    
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
        
class PartGlyph():
    
    def __init__(self, part):
        polydata = part.get_vtk()
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInput(polydata)
        actor = vtk.vtkActor()
        actor.GetProperty().SetColor(0,0,1)
        actor.SetMapper(mapper)
        self.actor = actor
        
        
class MeshGlyph():
    
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

class CellGlyph():
    
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
        
class VesselGlyph():
    
    def __init__(self, vessel_network):
        pointColors = vtk.vtkUnsignedCharArray()
        pointColors.SetNumberOfComponents(3)
        pointColors.SetName("PointColors") 
         
        converter = casie.rwc.vtk_converter.VesselConverter()
        vessel_vtk = converter.convert(vessel_network)
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