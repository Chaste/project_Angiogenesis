import logging
import os
import vtk
import wx
import numpy as np
from vmtk import vmtkscripts

import chaste.gui.properties
import chaste.plot.two.glyphs
import chaste.utility.rwc

if chaste.gui.properties._have_wx:
    import chaste.gui.panels.base
    
    class SurfaceFromImage2dPanel(chaste.gui.panels.base.Panel):
        
        ''' 
        Default panel
        
        Attributes
        ----------
        
        '''  
        
        def __init__(self, parent):
            
            ''' 
            Set up the panel, add the controls
            '''
            
            chaste.gui.panels.base.Panel.__init__(self, parent)
            self.name = "SurfaceFromImage2d"
            
        def add_controls(self):
            
            ''' 
            Add wx controls to the panel
            '''  
            
            self.select_input_file = wx.Button(parent = self, label="Load Surface File")
            self.file_label = wx.StaticText(parent = self, label="Input File: ")
            self.file_text = wx.StaticText(parent = self, label = "None")
            self.save_input_file = wx.Button(parent = self, label="Save Boundary File")
    
        def size_controls(self):
            
            ''' 
            Size the controls
            '''  
            
            centre_fmt = [0, wx.CENTER, 3]
            vbox = wx.BoxSizer(wx.VERTICAL)
            vbox.AddSpacer(10)
            vbox.Add(self.select_input_file, *centre_fmt)
            vbox.Add(self.file_label, *centre_fmt)
            vbox.Add(self.file_text, *centre_fmt)
            vbox.Add(self.save_input_file, *centre_fmt)
            vbox.AddSpacer(10)
    
            self.SetSizer(vbox)
            vbox.Fit(self) 
            
        def bind_events(self):
            
            ''' 
            Bind the events
            '''  
            
            self.select_input_file.Bind(wx.EVT_BUTTON, self.on_load_file)
            self.save_input_file.Bind(wx.EVT_BUTTON, self.on_save_file)
            
        def on_load_file(self, event = None):
            
            self.file_name = self.get_file_name()
            self.canvas = self.GetTopLevelParent().get_2d_canvas(show = True)
            self.canvas.add_tiff(self.file_name)
             
            logging.info("Loaded Image File: " + str(self.file_name))
            
            self.setup_tool()
            self.run_tool()
            
        def on_save_file(self, event = None):
            
            file_mod = os.path.splitext(self.file_name)[0] + "_surface.vtp"
            chaste.utility.rwc.write_vtk_surface(file_mod, self.surface)
            
            file_mod = os.path.splitext(self.file_name)[0] + "_boundaries.vtp"
            chaste.utility.rwc.write_vtk_surface(file_mod, self.boundaries)
            
            logging.info("Saved Surface File to : " + str(file_mod))
            
        def setup_tool(self):
            
            self.tool = SurfaceFromImage2d()
            image = chaste.utility.rwc.tiff_to_vti(self.file_name)
            self.tool.set_image(image)
            
        def run_tool(self):
            
            self.tool.update()
            self.surface = self.tool.get_output(surface = True)
            self.boundaries = self.tool.get_output(surface = False)
            
            glyph = chaste.visualization.two.glyphs.VtkLinesGlyph(self.boundaries)
            self.canvas.add_glyph(glyph, clear = False)
            
            logging.info("Extracted Surfaces")

class SurfaceFromImage2d():
    
    def __init__(self):
        
        self.surface = None
        self.boundaries = None
        self.image = None
        self.target_reduction = 0.9994
        self.num_subdivisions = 3
        self.resize_factors = None
        self.clipping_factor = 0.0
        
    def set_image(self, image):
        
        self.image = image
    
    def update(self):
        
        # Size the image if needed
        if self.resize_factors is not None:
            reduction_filter = vtk.vtkImageResize()
            reduction_filter.SetInput(self.image)
            reduction_filter.SetMagnificationFactors(self.resize_factors[0], self.resize_factors[1], self.resize_factors[2])
            reduction_filter.SetResizeMethodToMagnificationFactors()
            reduction_filter.Update()
        
        # Threshold
        threshold = vtk.vtkThreshold()
        if self.resize_factors is not None:
            threshold.SetInputConnection(reduction_filter.GetOutputPort())
        else:
            threshold.SetInput(self.image)
        threshold.ThresholdByUpper(1.0)
        threshold.Update()
        
        # Convert to polydata
        surface = vtk.vtkGeometryFilter()
        surface.SetInputConnection(threshold.GetOutputPort())
        surface.Update()
        
        # Triangulate
        triangle = vtk.vtkTriangleFilter()
        triangle.SetInputConnection(surface.GetOutputPort())
        triangle.Update()
    
        # Decimate
        decimate = vtk.vtkDecimatePro()
        decimate.SetInputConnection(triangle.GetOutputPort())
        decimate.SetTargetReduction(self.target_reduction)
        decimate.SetFeatureAngle(15.0)
        decimate.Update()
        
        # Do loop subdivision
        su = vtk.vtkLinearSubdivisionFilter()
        su.SetInputConnection(decimate.GetOutputPort())
        su.SetNumberOfSubdivisions(self.num_subdivisions)
        su.Update()
        
        # Clip the boundaries, recommended to ensure straight inlets and outlets
        if self.clipping_factor > 0.0:
            bounds = su.GetOutput().GetBounds()
            
            width = bounds[1] - bounds[0]
            height = bounds[3] - bounds[2]
            
            p1 = vtk.vtkPlane()
            p1.SetOrigin(bounds[0] + self.clipping_factor*width, 0.0, 0)
            p1.SetNormal(1, 0, 0)
            
            p2 = vtk.vtkPlane()
            p2.SetOrigin(0.0, bounds[2] + self.clipping_factor*height, 0)
            p2.SetNormal(0, 1, 0)
            
            p3 = vtk.vtkPlane()
            p3.SetOrigin(bounds[1] - self.clipping_factor*width, 0.0, 0)
            p3.SetNormal(-1, 0, 0)
        
            p4 = vtk.vtkPlane()
            p4.SetOrigin(0.0, bounds[3] - self.clipping_factor*height, 0)
            p4.SetNormal(0, -1, 0)        
        
            c1 = vtk.vtkClipPolyData()
            c1.SetInputConnection(su.GetOutputPort())
            c1.SetClipFunction(p1)
            c1.SetValue(0.0)
            
            c2 = vtk.vtkClipPolyData()
            c2.SetInputConnection(c1.GetOutputPort())
            c2.SetClipFunction(p2)
            c2.SetValue(0.0)
# 
            c3 = vtk.vtkClipPolyData()
            c3.SetInputConnection(c2.GetOutputPort())
            c3.SetClipFunction(p3)
            c3.SetValue(0.0)
             
            c4 = vtk.vtkClipPolyData()
            c4.SetInputConnection(c3.GetOutputPort())
            c4.SetClipFunction(p4)
            c4.SetValue(0.0)
            
            su = vtk.vtkGeometryFilter()
            su.SetInputConnection(c4.GetOutputPort())
            su.Update()
        
        feature_edges = vtk.vtkFeatureEdges()
        feature_edges.SetInputConnection(su.GetOutputPort())  
        feature_edges.Update()   
        
        clean = vtk.vtkCleanPolyData()
        clean.SetInputConnection(feature_edges.GetOutputPort())
        clean.Update()
         
        triangle2 = vtk.vtkTriangleFilter()
        triangle2.SetInputConnection(clean.GetOutputPort())
        triangle2.Update()
        
        self.surface = su.GetOutput()
        self.boundaries = triangle2.GetOutput()
        self.boundaries = self.boundary_to_polylines()
        
        smoother = vmtkscripts.vmtkCenterlineResampling()
        smoother.Centerlines = self.boundaries
        smoother.length = 3.0
        smoother.Execute()
        self.boundaries = smoother.Centerlines
        
    def boundary_to_polylines(self):
        
        numPoints = self.boundaries.GetNumberOfPoints()    
        vtkpoints = self.boundaries.GetPoints()  
        connectivity = []
        
        points = [] 
        for i in range(numPoints):
            points.append(vtkpoints.GetPoint(i))
            connectivity.append([])
            
        numCells = self.boundaries.GetNumberOfLines()  
        cellArray = self.boundaries.GetLines()
        cellArray.InitTraversal()
        segList = vtk.vtkIdList()
            
        edges = []
        for i in range(numCells): 
            cellArray.GetNextCell(segList)
            point_indices = []
            for j in range(0, segList.GetNumberOfIds()):
                seg_id = segList.GetId(j)
                point_indices.append(int(seg_id))
            edges.append((point_indices[0], point_indices[1]))
            connectivity[point_indices[0]].append(i)
            connectivity[point_indices[1]].append(i)
            
        regions = []
        point_visited = np.zeros(len(points))
        
        for idx, eachPoint in enumerate(points):
            if point_visited[idx] == 0:
                point_visited[idx] = 1
                region_points = [idx]
                current_point_id = idx
                previous_point_id = idx
                found_visited = False
                while not found_visited:
                    edge_1 = edges[connectivity[current_point_id][0]]
                    edge_2 = edges[connectivity[current_point_id][1]]
                    if edge_1[0] == current_point_id:
                        opp1 = edge_1[1]
                    else:
                        opp1 = edge_1[0]
                    if opp1 != previous_point_id:
                        next_edge = edge_1
                    else:
                        next_edge = edge_2
                              
                    if next_edge[0]== current_point_id:
                        next_point_id = next_edge[1]  
                    else:
                        next_point_id = next_edge[0]     
                              
                    if point_visited[next_point_id] == 1:
                        region_points.append(next_point_id)
                        break
                      
                    point_visited[next_point_id] = 1
                    region_points.append(next_point_id)
                    previous_point_id = current_point_id
                    current_point_id = next_point_id
                          
                if len(region_points) > 1:
                    regions.append(region_points)
                        
        new_points = vtk.vtkPoints()
        new_points.SetNumberOfPoints(len(points))
        for idx, eachPoint in enumerate(points):
            new_points.SetPoint(idx, points[idx][0], points[idx][1], 0.0)
            
        lines = vtk.vtkCellArray()
        for eachRegion in regions:
            lines.InsertNextCell(len(eachRegion))
            for eachPoint in eachRegion:
                lines.InsertCellPoint(eachPoint)
                
        polygon = vtk.vtkPolyData()
        polygon.SetPoints(new_points)
        polygon.SetLines(lines)
        
        return polygon  

    def get_output(self, surface = True):
        if surface:
            return self.surface
        else:
            return self.boundaries