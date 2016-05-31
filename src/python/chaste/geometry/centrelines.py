import os
import logging
import numpy as np
import vtk
from vmtk import vmtkscripts
from scipy.spatial import Voronoi
from scipy.sparse import csr_matrix
import scipy.sparse.csgraph
import wx

import casie.gui.properties
import casie.geometry.converters
import casie.plot.two.glyphs
import casie.utility.rwc
import casie.population.vessel

if casie.gui.properties._have_wx:
    import casie.gui.panels.base
    
    class Centrelines2dPanel(casie.gui.panels.base.Panel):
        
        ''' 
        Default panel
        
        Attributes
        ----------
        
        '''  
        
        def __init__(self, parent):
            
            ''' 
            Set up the panel, add the controls
            '''
            
            casie.gui.panels.base.Panel.__init__(self, parent)
            self.name = "Centrelines2d"
            
        def add_controls(self):
            
            ''' 
            Add wx controls to the panel
            '''  
            
            self.select_input_file = wx.Button(parent = self, label="Load Surface File")
            self.file_label = wx.StaticText(parent = self, label="Input File: ")
            self.file_text = wx.StaticText(parent = self, label = "None")
            self.save_input_file = wx.Button(parent = self, label="Save File")
    
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
            self.surface = casie.utility.rwc.read_vtk_surface(self.file_name)
     
            glyph = casie.plot.two.glyphs.VtkLinesGlyph(self.surface)
            
            self.canvas = self.GetTopLevelParent().get_2d_canvas(show = False)
            self.canvas.add_glyph(glyph, True)
             
            logging.info("Loaded Boundary File: " + str(self.file_name))
            
            self.setup_tool()
            self.run_tool()
            
        def on_save_file(self, event = None):
            
            file_mod = os.path.splitext(self.file_name)[0] + "_centres.vtp"
            self.network.write(str(file_mod))
            
            logging.info("Saved Centreline File to : " + str(file_mod))
            
        def setup_tool(self):
            
            self.tool = Centrelines2d()
            
        def run_tool(self):
            
            self.tool.set_surface(self.surface)
            self.tool.update()
            self.network = self.tool.get_output()
            
            glyph = casie.plot.two.glyphs.VesselNetworkGlyph(self.network)
            self.canvas.add_glyph(glyph, clear = False)
            
            logging.info("Extracted Centreline")

class Centrelines2d():
    
    def __init__(self):
        self.surface = None
        self.network = None
        self.start_points= None
        
    def set_surface(self, surface):
        
        self.surface = surface

    def update(self):
        
        # Convert from vtk format to tri plc format
        vtk_to_tri = casie.geometry.converters.VtkToTri()
        points, edges = vtk_to_tri.generate(self.surface)
        
        # Get the voronoi diagram of the input points
        input_points = []
        for eachEdge in edges:
            input_points.append([eachEdge[0], eachEdge[1]])
        vor = Voronoi(points)
        verts = vor.vertices
        ridges = vor.ridge_vertices
        
        # Get the finite ridges
        finite_ridges = []
        for eachRidge in ridges:
            if eachRidge[0] >= 0 and eachRidge[1] >= 0:
                finite_ridges.append(eachRidge)
        finite_ridges
        
        # Only keep ridges which do not intersect the original polygon
        internal_ridges = []
        for eachRidge in finite_ridges:
            
            crosses = False
            for eachEdge in edges:
                p1 = (points[eachEdge[0]][0], points[eachEdge[0]][1], 0.0)
                p2 = (points[eachEdge[1]][0], points[eachEdge[1]][1], 0.0)
                x1 = (verts[eachRidge[0]][0], verts[eachRidge[0]][1], 0.0)
                x2 = (verts[eachRidge[1]][0], verts[eachRidge[1]][1], 0.0)
    
                para1 = vtk.mutable(0)
                para2 = vtk.mutable(0)
                intersect = vtk.vtkLine.Intersection(p1, p2, x1, x2, para1, para2)
                if intersect != 0:
                    crosses = True
                    break
            if not crosses:
                internal_ridges.append(eachRidge)
                
        # Find the largest connected path
        G_dense = np.zeros((len(verts), len(verts)))
        G_verts = []
        for idx in range(len(verts)):
            G_verts.append([])
        
        for eachRidge in internal_ridges:
            G_verts[eachRidge[0]].append(eachRidge[1])
            G_verts[eachRidge[1]].append(eachRidge[0])
            G_dense[eachRidge[0], eachRidge[1]] = 1.0
            G_dense[eachRidge[1], eachRidge[0]] = 1.0
            
        G_sparse = csr_matrix(G_dense)
        num, labels = scipy.sparse.csgraph.connected_components(G_sparse, False)
        label_count = [np.sum(labels == i) for i in range(num)]
        max_label = label_count.index(max(label_count))
        sparse_edge_labels =  np.where(labels == max_label)
        
        sparse_edges = []
        for eachLabel in sparse_edge_labels[0]:
            for eachNeighbour in G_verts[eachLabel]:
                sparse_edges.append([eachLabel, eachNeighbour])
                
        # Use VTK to remove duplicate points
        tri_to_vtk = casie.geometry.converters.TriToVtk([verts, sparse_edges])
        polyData = tri_to_vtk.generate()
        
        clean = vtk.vtkCleanPolyData()
        clean.SetInput(polyData)
        clean.Update()
        
        casie.utility.rwc.write_vtk_surface("/home/grogan/voronoi.vtp", clean.GetOutput())  
        
        # Remove duplicate edges
        centre = clean.GetOutput()
        vtk_to_tri = casie.geometry.converters.VtkToTri()
        points, edges = vtk_to_tri.generate(centre)
        unique_edges = []

        for idx, eachEdge in enumerate(edges):
            duplicate = False
            for eachInnerEdge in unique_edges:
                if eachInnerEdge[0] == eachEdge[0] and eachInnerEdge[1] == eachEdge[1]:
                    duplicate = True
                    break
                if eachInnerEdge[1] == eachEdge[0] and eachInnerEdge[0] == eachEdge[1]:
                    duplicate = True
                    break 
            if not duplicate:
                unique_edges.append(eachEdge)
        new_vtk_lines = vtk.vtkCellArray()
        for eachEdge in unique_edges:
            line = vtk.vtkLine()
            line.GetPointIds().InsertId(0, eachEdge[0])
            line.GetPointIds().InsertId(1, eachEdge[1])
            new_vtk_lines.InsertNextCell(line)
        
        tidy_edges = vtk.vtkPolyData() 
        tidy_edges.SetPoints(centre.GetPoints())
        tidy_edges.SetLines(new_vtk_lines)
        
        clean = vtk.vtkCleanPolyData()
        clean.SetInput(tidy_edges)
        clean.Update()
        tidy_edges = clean.GetOutput()
        casie.utility.rwc.write_vtk_surface("/home/grogan/voronoi_tidy.vtp", tidy_edges)  
        
        vtk_numPoints = tidy_edges.GetNumberOfPoints()    
        vtk_points = tidy_edges.GetPoints() 
        cellArray = tidy_edges.GetLines()
        locator = vtk.vtkKdTreePointLocator()
        locator.SetDataSet(tidy_edges)
        
        for idx, eachLoc in enumerate(self.start_points):
            probe_loc = np.array((eachLoc[0], eachLoc[1], 0.0))
            closest_id = locator.FindClosestPoint(probe_loc)
            vtk_points.InsertNextPoint(probe_loc)
            
            cellArray.InsertNextCell(2)
            cellArray.InsertCellPoint(closest_id)
            cellArray.InsertCellPoint(vtk_numPoints + idx)
        polygon = vtk.vtkPolyData()
        polygon.SetPoints(vtk_points)
        polygon.SetLines(cellArray)
        
        clean = vtk.vtkCleanPolyData()
        clean.SetInput(polygon)
        clean.Update()
        polygon = clean.GetOutput()
        casie.utility.rwc.write_vtk_surface("/home/grogan/voronoi_extd.vtp", polygon)  
        
        pruned = self.prune(polygon)
        pruned = self.prune(pruned)
        casie.utility.rwc.write_vtk_surface("/home/grogan/voronoi_pruned.vtp", pruned) 
        
        polybound = casie.geometry.converters.vtk_network_lines_to_polylines(pruned)
        
        smoother = vmtkscripts.vmtkCenterlineResampling()
        smoother.Centerlines = polybound
        smoother.length = 15.0
        smoother.Execute()
        self.network = smoother.Centerlines
        
        boundary_point_label = vtk.vtkFloatArray()
        boundary_point_label.SetNumberOfComponents(1)
        boundary_point_label.SetName("Radius")    
        
        locator = vtk.vtkKdTreePointLocator()
        locator.SetDataSet(self.surface) 
        vtk_numPoints = self.network.GetNumberOfPoints()    
        vtk_points = self.network.GetPoints() 
        surf_points = self.surface.GetPoints() 
        
        for idx in range(vtk_numPoints):
            loc = vtk_points.GetPoint(idx)
            closest = locator.FindClosestPoint(loc)
            surf_closest = surf_points.GetPoint(closest)
            radius = np.linalg.norm(np.array(loc) - np.array(surf_closest))/2.0
            boundary_point_label.InsertNextTupleValue((float(radius),))  
        
        self.network.GetPointData().AddArray(boundary_point_label)  

    def prune(self, surface):
        
        numPoints = surface.GetNumberOfPoints()    
        vtkpoints = surface.GetPoints()  
        connectivity = []
        points = [] 
        for i in range(numPoints):
            points.append(vtkpoints.GetPoint(i))
            connectivity.append([])
            
        numCells = surface.GetNumberOfLines()  
        cellArray = surface.GetLines()
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
            
        c1_points = np.zeros(numPoints)
        locator = vtk.vtkKdTreePointLocator()
        locator.SetDataSet(surface)
        near_seeds = np.zeros(numPoints)
        for idx, eachLoc in enumerate(self.start_points):
            probe_loc = np.array((eachLoc[0], eachLoc[1], 0.0))
            closest_id = locator.FindClosestPoint(probe_loc)
            near_seeds[closest_id] = 1.0
        
        for idx, eachPoint in enumerate(points):
            num_segs = len(connectivity[idx])
            if num_segs == 1 and near_seeds[idx]==0:
                c1_points[idx] = 1
                current_point_id = idx
                previous_point_id = idx
                found_branch = False
                while not found_branch:
                    current_num_segs = len(connectivity[current_point_id])
                    if current_num_segs == 2:
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
                    else:
                        next_edge = edges[connectivity[idx][0]]
                              
                    if next_edge[0]== current_point_id:
                        next_point_id = next_edge[1]  
                    else:
                        next_point_id = next_edge[0]     
                              
                    if len(connectivity[next_point_id]) == 2:
                        c1_points[next_point_id] = 1
                        previous_point_id = current_point_id
                        current_point_id = next_point_id
                    else:
                        found_branch = True
                        break
            
        boundary_point_label = vtk.vtkFloatArray()
        boundary_point_label.SetNumberOfComponents(1)
        boundary_point_label.SetName("PointRemoveLabel")     
        
        for idx, eachPoint in enumerate(c1_points):
            boundary_point_label.InsertNextTupleValue((float(c1_points[idx]),))  
            
        surface.GetPointData().AddArray(boundary_point_label)  
        
        threshold = vtk.vtkThreshold()
        threshold.SetInput(surface)
        threshold.SetInputArrayToProcess(0, 0, 0, "vtkDataObject::FIELD_ASSOCIATION_POINTS", "PointRemoveLabel")
        threshold.ThresholdBetween(0, 0)
        threshold.Update()
            
        polysurface = vtk.vtkGeometryFilter()
        polysurface.SetInputConnection(threshold.GetOutputPort())
        polysurface.Update()
        
        return polysurface.GetOutput()
        
    def get_output(self):
        return self.network
        