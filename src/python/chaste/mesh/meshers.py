import logging
import os
import collections
import vtk
import wx
import dolfin as df
import numpy as np
from meshpy.triangle import MeshInfo, build

import chaste.utility.rwc
import chaste.geometry.labelling
import chaste.geometry.converters.other
import chaste.plot.two.glyphs
import chaste.gui.properties
import chaste.population.vessel
import chaste.mesh.centrelines
   
if chaste.gui.properties._have_wx:
    import chaste.gui.panels.base
 
    class Mesh2dPanel(chaste.gui.panels.base.Panel):
     
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
            self.name = "Mesh2d"
            self.surface = None
            self.domain = None
            self.centres = None
            self.region1_points = None
            self.region2_points = None
            self.holes = None
            self.picking_r1 = True
            self.picking_r2 = False
            self.picking_holes = False
            self.canvas_set = False
            self.boundary_edges = None
            self.domain_edges = None
            self.lines = None
            self.network = None
            
        def add_controls(self):
            
            ''' 
            Add wx controls to the panel
            '''  
            
            self.select_input_file = wx.Button(parent = self, label="Load Network Surface File")
            self.select_input_domain_file = wx.Button(parent = self, label="Load Domain File")
            self.select_input_centreline_file = wx.Button(parent = self, label="Load Network Centreline File")
            self.region1_mesh_size = wx.StaticText(parent = self, label = "Region 1 Mesh Size")
            self.r1_mesh = wx.TextCtrl(self, value="5.0")
            self.region2_mesh_size = wx.StaticText(parent = self, label = "Region 2 Mesh Size")
            self.r2_mesh = wx.TextCtrl(self, value="15.0")
            self.holes = wx.Button(parent = self, label="Pick Holes")
            self.r1_point_pick = wx.Button(parent = self, label="Pick Region 1 Points")
            self.r2_point_pick = wx.Button(parent = self, label="Pick Region 2 Points")
            self.generate = wx.Button(parent = self, label="Generate")
            self.save_input_file = wx.Button(parent = self, label="Save Mesh File")

        def size_controls(self):
            
            ''' 
            Size the controls
            '''  
            
            centre_fmt = [0, wx.CENTER, 3]
            vbox = wx.BoxSizer(wx.VERTICAL)
            vbox.AddSpacer(10)
            vbox.Add(self.select_input_file, *centre_fmt)
            vbox.Add(self.select_input_domain_file, *centre_fmt)
            vbox.Add(self.select_input_centreline_file, *centre_fmt)
            vbox.Add(self.region1_mesh_size, *centre_fmt)
            vbox.Add(self.r1_mesh, *centre_fmt)
            vbox.Add(self.region2_mesh_size, *centre_fmt)
            vbox.Add(self.r2_mesh, *centre_fmt)
            vbox.Add(self.holes, *centre_fmt)
            vbox.Add(self.r1_point_pick, *centre_fmt)
            vbox.Add(self.r2_point_pick, *centre_fmt)
            vbox.Add(self.generate, *centre_fmt)
            vbox.Add(self.save_input_file, *centre_fmt)
            vbox.AddSpacer(10)
        
            self.SetSizer(vbox)
            vbox.Fit(self) 
            
        def bind_events(self):
            
            ''' 
            Bind the events
            '''  
            
            self.select_input_file.Bind(wx.EVT_BUTTON, self.on_load_file)
            self.select_input_domain_file.Bind(wx.EVT_BUTTON, self.on_load_domain_file)
            self.select_input_centreline_file.Bind(wx.EVT_BUTTON, self.on_load_centreline_file)
            self.holes.Bind(wx.EVT_BUTTON, self.on_pick_holes)
            self.r1_point_pick.Bind(wx.EVT_BUTTON, self.on_pick_r1)
            self.r2_point_pick.Bind(wx.EVT_BUTTON, self.on_pick_r2)
            self.generate.Bind(wx.EVT_BUTTON, self.on_generate)
            self.save_input_file.Bind(wx.EVT_BUTTON, self.on_save_file)
            
        def on_pick_holes(self, event = None):
            self.holes = []
            self.picking_r1 = False
            self.picking_r2 = False
            self.picking_holes = True
            
        def on_pick_r2(self, event = None):
            self.region2_points = []
            self.picking_r1 = False
            self.picking_r2 = True
            self.picking_holes = False

        def on_pick_r1(self, event = None):
            self.region1_points = []
            self.picking_r1 = True
            self.picking_r2 = False
            self.picking_holes = False            
            
        def on_load_file(self, event = None):
            
            self.file_name = self.get_file_name()
            self.surface = chaste.utility.rwc.read_vtk_surface(self.file_name, True, True)
            
            bound_extractor = chaste.geometry.labelling.BoundaryExtractor2d()
            bound_extractor.labels = [1,2]
            bound_extractor.surface = self.surface
            bound_extractor.update()
            self.boundary_edges = bound_extractor.get_output()
            
            glyph = chaste.plot.two.glyphs.VtkLinesGlyph(self.surface, color = "red")
            self.canvas = self.GetTopLevelParent().get_2d_canvas(show = False)
            self.canvas.add_glyph(glyph, self.domain is None)
             
            if not self.canvas_set:
                self.setup_canvas()
            logging.info("Loaded Boundary File: " + str(self.file_name))
            
        def on_load_domain_file(self, event = None):
            
            self.file_name = self.get_file_name()
            self.domain = chaste.utility.rwc.read_vtk_surface(self.file_name)
            self.lines = chaste.geometry.other.vtkpolygon_to_lines(self.domain)
        
            glyph = chaste.plot.two.glyphs.VtkLinesGlyph(self.lines, color="blue")
            self.canvas = self.GetTopLevelParent().get_2d_canvas(show = False)
            self.canvas.add_glyph(glyph, True)
            
            if not self.canvas_set:
                self.setup_canvas()
                
            conv = chaste.geometry.other.VtkToTri()
            points, edges = conv.generate(self.lines)
            myedges = []
            for eachEdge in edges:
                myedges.append([points[eachEdge[0]], points[eachEdge[1]]])
            self.domain_edges = myedges
             
            logging.info("Loaded Domain File: " + str(self.file_name))
            
        def on_load_centreline_file(self, event = None):
            
            self.file_name = self.get_file_name()
            self.centres = chaste.utility.rwc.read_vtk_surface(self.file_name, True, True)
            self.network = self.centres
            glyph = chaste.plot.two.glyphs.VtkLinesGlyph(self.centres, color="green")
            self.canvas = self.GetTopLevelParent().get_2d_canvas(show = False)
            self.canvas.add_glyph(glyph, self.domain is None)
            
            if not self.canvas_set:
                self.setup_canvas()
             
            logging.info("Loaded Domain File: " + str(self.file_name))
            
        def setup_canvas(self):
            
            self.canvas.canvas.mpl_connect('button_press_event', self.on_press)
            
        def on_press(self, event):
            
            self.canvas.temp_sketch.new_point([event.xdata, event.ydata])
            if self.picking_r1:
                self.region1_points.append([event.xdata, event.ydata])
            elif self.picking_r2:
                self.region2_points.append([event.xdata, event.ydata])
            else:
                self.holes.append([event.xdata, event.ydata])
            self.canvas.re_draw()   
            
        def on_generate(self, event = None):
            
            logging.info("Generating Mesh")
            self.setup_tool()
            self.run_tool()
            logging.info("Mesh Generated")
            
        def on_save_file(self, event = None):
            
            converter = DolfinConverter2d()
            converter.boundary_edges = self.boundary_edges
            converter.domain_edges = self.domain_edges
            converter.network = self.network

            main_mesh, r1_mesh, r2_mesh = converter.update(self.mesh, regions_set=self.tool.region1_points is not None)
            
            if converter.new_network is not None:
                file_mod = os.path.splitext(self.file_name)[0] + "_new_network.vtp"
                chaste.utility.rwc.write_vtk_surface(file_mod, converter.new_network)
            
            file_mod = os.path.splitext(self.file_name)[0] + "_mesh.pvd"
            path = str(file_mod)
            outfile = df.File(path)
            outfile << main_mesh[0]
            
            file_mod = os.path.splitext(self.file_name)[0] + "_mesh.xml"
            path = str(file_mod)
            outfile = df.File(path)
            outfile << main_mesh[0]
            
            file_mod = os.path.splitext(self.file_name)[0] + "_mesh_domains.pvd"
            path = str(file_mod)
            outfile = df.File(path)
            outfile << main_mesh[1]
            
            file_mod = os.path.splitext(self.file_name)[0] + "_mesh_domains.xml"
            path = str(file_mod)
            outfile = df.File(path)
            outfile << main_mesh[1]
            
            file_mod = os.path.splitext(self.file_name)[0] + "_mesh_boundaries.pvd"
            path = str(file_mod)
            outfile = df.File(path)
            outfile << main_mesh[2]    
            
            file_mod = os.path.splitext(self.file_name)[0] + "_mesh_boundaries.xml"
            path = str(file_mod)
            outfile = df.File(path)
            outfile << main_mesh[2]    
            
            if r1_mesh is not None:
                file_mod = os.path.splitext(self.file_name)[0] + "_mesh_r1.pvd"
                path = str(file_mod)
                outfile = df.File(path)
                outfile << r1_mesh[0]
                
                file_mod = os.path.splitext(self.file_name)[0] + "_mesh_r1.xml"
                path = str(file_mod)
                outfile = df.File(path)
                outfile << r1_mesh[0]
                
                file_mod = os.path.splitext(self.file_name)[0] + "_mesh_r1_domains.pvd"
                path = str(file_mod)
                outfile = df.File(path)
                outfile << r1_mesh[1]
                
                file_mod = os.path.splitext(self.file_name)[0] + "_mesh_r1_domains.xml"
                path = str(file_mod)
                outfile = df.File(path)
                outfile << r1_mesh[1]
                
                file_mod = os.path.splitext(self.file_name)[0] + "_mesh_r1_boundaries.pvd"
                path = str(file_mod)
                outfile = df.File(path)
                outfile << r1_mesh[2]   
                
                file_mod = os.path.splitext(self.file_name)[0] + "_mesh_r1_boundaries.xml"
                path = str(file_mod)
                outfile = df.File(path)
                outfile << r1_mesh[2]  
                
            if r2_mesh is not None: 
            
                file_mod = os.path.splitext(self.file_name)[0] + "_mesh_r2.pvd"
                path = str(file_mod)
                outfile = df.File(path)
                outfile << r2_mesh[0]  
                
                file_mod = os.path.splitext(self.file_name)[0] + "_mesh_r2.xml"
                path = str(file_mod)
                outfile = df.File(path)
                outfile << r2_mesh[0]  
                
                file_mod = os.path.splitext(self.file_name)[0] + "_mesh_r2_domains.pvd"
                path = str(file_mod)
                outfile = df.File(path)
                outfile << r2_mesh[1]
                
                file_mod = os.path.splitext(self.file_name)[0] + "_mesh_r2_domains.xml"
                path = str(file_mod)
                outfile = df.File(path)
                outfile << r2_mesh[1]
                
                file_mod = os.path.splitext(self.file_name)[0] + "_mesh_r2_boundaries.pvd"
                path = str(file_mod)
                outfile = df.File(path)
                outfile << r2_mesh[2]    
                
                file_mod = os.path.splitext(self.file_name)[0] + "_mesh_r2_boundaries.xml"
                path = str(file_mod)
                outfile = df.File(path)
                outfile << r2_mesh[2]              
            
            logging.info("Saved Mesh File to : " + str(file_mod))
            
        def setup_tool(self):
            
            self.tool = Mesher2d()
            if self.surface is not None:
                self.tool.set_vtk_surface(self.surface)
            if self.domain is not None:
                self.tool.domain = self.lines
            if self.centres is not None:
                self.tool.centres = self.centres
                
            self.tool.holes = self.holes
            self.tool.region1_points = self.region1_points
            self.tool.region2_points = self.region2_points
            self.tool.region1_mesh_size = float(self.r1_mesh.GetValue())
            self.tool.region2_mesh_size = float(self.r2_mesh.GetValue())
            
        def run_tool(self):
            
            self.tool.update()
            self.mesh = self.tool.get_output()
            
            glyph = chaste.plot.two.glyphs.MeshGlyph([self.mesh.points, self.mesh.elements], False)
            self.canvas.add_glyph(glyph, clear = False)
            
            logging.info("Generate Mesh")

class Mesher2d():
    
    def __init__(self):
        
        self.mesh = None
        self.region_markers = None
        self.domain = None
        self.network = None
        self.centres = None
        self.region_markers = None
        self.surface = True
        self.region1_mesh_size = 10.0
        self.region2_mesh_size = 10.0
        self.region1_points = None
        self.region2_points = None
        self.holes = None
        self.vtk_surface = None
        
    def set_vtk_surface(self, surface):
        
        self.vtk_surface = surface
        
    def vtk_to_points(self, polydata, points, edges):
        
        num_domain_points = len(points)
        
        vtk_numPoints = polydata.GetNumberOfPoints()    
        vtk_points = polydata.GetPoints()   
        for i in range(vtk_numPoints):
            points.append([vtk_points.GetPoint(i)[0], vtk_points.GetPoint(i)[1]])
               
        numCells = polydata.GetNumberOfLines()  
        cellArray = polydata.GetLines()
        cellArray.InitTraversal()
        segList = vtk.vtkIdList()
        
        for i in range(numCells): 
            cellArray.GetNextCell(segList)
            point_indices = []
            for j in range(0, segList.GetNumberOfIds()):
                seg_id = segList.GetId(j)
                point_indices.append(num_domain_points + int(seg_id))
            edges.append(point_indices)
        
        return points, edges
    
    def update(self):
        
        # If no vessel surface is required can directly use points and edges
        points = []
        edges = []
        
        if self.domain is not None:
            points, edges = self.vtk_to_points(self.domain, points, edges)
            
        if self.centres is not None:
            points, edges = self.vtk_to_points(self.centres, points, edges)        
            
#         if self.network is not None:
#             if not self.surface:
#                 num_points = len(points)
#                 nodes = self.network.nodes
#     
#                 for idx, eachNode in enumerate(nodes):
#                     points.append([eachNode.GetLocationVector()[0],eachNode.GetLocationVector()[1]] )
#                     eachNode.id = idx + num_points
#                     
#                 vessels = self.network.vessels
#                 
#                 for eachVessel in vessels:
#                     edges.append((eachVessel.start_node.id, eachVessel.end_node.id))
#             
#             else:
#                 converter = chaste.geometry.converters.other.NetworkToGeometry(self.network, dimension = 2)
#                 network_edges = converter.generate()
#                 
#                 points, edges = self.vtk_to_points(network_edges, points, edges)
                    
        if self.vtk_surface is not None and self.centres is None:
            points, edges = self.vtk_to_points(self.vtk_surface, points, edges)
            
        # Do the meshing with triangle
        mesh_info = MeshInfo()
        mesh_info.set_points(points)
        mesh_info.set_facets(edges)
        
        if self.region1_points is not None:
            
            total_regions = len(self.region1_points)
            if self.region2_points is not None:
                total_regions += len(self.region2_points)
            
            mesh_info.regions.resize(total_regions)
            for idx, eachPoint in enumerate(self.region1_points):
                mesh_info.regions[idx] = [eachPoint[0],eachPoint[1],1,self.region1_mesh_size]
            if self.region2_points is not None:
                for idx, eachPoint in enumerate(self.region2_points):
                    mesh_info.regions[idx + len(self.region1_points)] = [eachPoint[0],eachPoint[1],2,self.region2_mesh_size]
            if self.holes is not None:
                mesh_info.holes.resize(len(self.holes))
                for idx, eachHole in enumerate(self.holes):
                    mesh_info.holes[idx] = [eachHole[0], eachHole[1]]
                            
            self.mesh = build(mesh_info, volume_constraints=True, attributes=True, generate_faces=True)
        else:
            
            if self.holes is not None:
                mesh_info.holes.resize(len(self.holes))
                for idx, eachHole in enumerate(self.holes):
                    mesh_info.holes[idx] = [eachHole[0], eachHole[1]]
            self.mesh = build(mesh_info, max_volume = self.region1_mesh_size)
            
    def get_output(self):
        
        return self.mesh
    
class DolfinConverter2d():
    
    def __init__(self):
        
        self.boundary_edges = None
        self.domain_edges = None
        self.network = None
        self.surface = None
        self.new_network = None
        
    def update(self, mesh, regions_set = False):
        editor = df.MeshEditor()
        dolfin_mesh = df.Mesh()
        
        editor.open(dolfin_mesh, 2, 2)
        editor.init_vertices(len(mesh.points))
        editor.init_cells(len(mesh.elements))
        
        print len(mesh.points), len(mesh.elements)
        
        element_labels = []
        for i, p in enumerate(mesh.points):
            editor.add_vertex(i, np.array(p))
        
        for i, t in enumerate(mesh.elements):
            editor.add_cell(i, np.array(t, dtype=np.uintp))
            if regions_set:
                element_labels.append(mesh.element_attributes[i])
        editor.close()
        dolfin_mesh.init()
        
                # Mark the domains
        mesh_func = df.MeshFunction("size_t", dolfin_mesh, 2)
        mesh_func.set_all(0)
        
        if len(element_labels)>0:
            num_cells = len(dolfin_mesh.cells())
            for idx in range(num_cells):
                mesh_func[idx] = int(element_labels[idx])
    
        # Set up facet markers
        boundaries = df.FacetFunction("size_t", dolfin_mesh)
        boundaries.set_all(0)
        
        if self.network is not None:
            self.update_centrelines(dolfin_mesh, boundaries)
            main_mesh_data = [dolfin_mesh, mesh_func, boundaries]
            
            start_points = []
            for eachRegion in self.boundary_edges:
                for eachEdge in eachRegion:
                    start_points.append((np.array(eachEdge[0]) + np.array(eachEdge[1]))/2.0)
            
            centre_tool = chaste.mesh.centrelines.Extract2d(dolfin_mesh, boundaries, 3)
            centre_tool.seed_points = start_points
            centre_mesh, new_network = centre_tool.generate()
            self.new_network = new_network
            
            r1mesh_func = df.MeshFunction("size_t", centre_mesh, 2)
            r1mesh_func.set_all(0)
            
            r1boundaries = df.FacetFunction("size_t", centre_mesh)
            r1boundaries.set_all(0)
            
            r1_mesh_data = [centre_mesh, r1mesh_func, r1boundaries]
            
            return main_mesh_data, r1_mesh_data, None
        
        facet_count = 0
        two_regions_found = False
        for eachFacet in df.facets(dolfin_mesh):
            if eachFacet.num_entities(2) == 2:
                if mesh_func[int(eachFacet.entities(2)[0])] != mesh_func[int(eachFacet.entities(2)[1])]: 
                    boundaries[facet_count] = 3
                    two_regions_found = True
                elif mesh_func[int(eachFacet.entities(2)[0])] == 1:
                    boundaries[facet_count] = 4
                elif mesh_func[int(eachFacet.entities(2)[0])] == 2: 
                    boundaries[facet_count] = 5
            facet_count += 1
         
        if two_regions_found: 
            r1mesh = df.SubMesh(dolfin_mesh, mesh_func, 1)
            r1mesh_func = df.MeshFunction("size_t", r1mesh, 2)
            r1mesh_func.set_all(0)
            
            r1boundaries = df.FacetFunction("size_t", r1mesh)
            r1boundaries.set_all(0)
            
            r2mesh = df.SubMesh(dolfin_mesh, mesh_func, 2)
            r2mesh_func = df.MeshFunction("size_t", r2mesh, 2)
            r2mesh_func.set_all(0) 
            
            r2boundaries = df.FacetFunction("size_t", r2mesh)
            r2boundaries.set_all(0)
        
        if self.boundary_edges is not None:
            line_marker = LineBoundary()
            line_marker.set_edges(self.boundary_edges[0])
            line_marker.mark(boundaries, 1)
            
            line_marker.set_edges(self.boundary_edges[1])    
            line_marker.mark(boundaries, 2)
            
            line_marker.set_edges(self.boundary_edges[0])
        
            if two_regions_found:
                default_boundary = DefaultBoundary()
                default_boundary.mark(r1boundaries, 3)   
                line_marker.mark(r1boundaries, 1)
                line_marker.set_edges(self.boundary_edges[1])    
                line_marker.mark(boundaries, 2)
                line_marker.mark(r1boundaries, 2)
    
                default_boundary = DefaultBoundary()
                default_boundary.mark(r2boundaries, 3)    
                
                if self.domain_edges is not None:
                    line_marker.set_edges(self.domain_edges)
                    line_marker.mark(r2boundaries, 0)
        
        main_mesh_data = [dolfin_mesh, mesh_func, boundaries]
        if two_regions_found:
            r1_mesh_data = [r1mesh, r1mesh_func, r1boundaries]
            r2_mesh_data = [r2mesh, r2mesh_func, r2boundaries]
        else:
            r1_mesh_data = None
            r2_mesh_data = None        
        
        return main_mesh_data, r1_mesh_data, r2_mesh_data
    
    def update_centrelines(self, mesh, boundaries):
        
        centres = OnVesselCentre()
        centres.set_network(self.network)
        centres.mark(boundaries, 3)
    
class LineBoundary(df.SubDomain):
    
    def set_edges(self, edges, tol = 1.e-3):
        self.edges = edges
        self.tol = tol
    
    def inside(self, x, on_boundary):
        
        inside = False
        
        if len(x) == 2:
            position = np.array((x[0], x[1], 0.0))
        else:
            position = np.array(x)
            
        if on_boundary:
            for eachEdge in self.edges:
                if len(eachEdge[0]) == 2:
                    e1loc = np.array((eachEdge[0][0], eachEdge[0][1], 0.0))
                    e2loc = np.array((eachEdge[1][0], eachEdge[1][1], 0.0))
                else:
                    e1loc = np.array(eachEdge[0])
                    e2loc = np.array(eachEdge[1])
                if vtk.vtkLine.DistanceToLine(position, e1loc, e2loc) <= self.tol:
                    dp1 = np.linalg.norm(e1loc - position)
                    dp2 = np.linalg.norm(e2loc - position)
                    dpLine = np.linalg.norm(e1loc - e2loc)
                    
                    if dp1 + dp2 <= dpLine + self.tol:
                        inside = True
                        break
        return inside
    
class DefaultBoundary(df.SubDomain):
    
    def inside(self, x, on_boundary):
        return on_boundary
    
class OnVesselCentre(df.SubDomain):
     
    def set_network(self, network, tol = 1.e-3):
        self.network = network
        self.tol = tol
        vtk_to_tri = chaste.geometry.other.VtkToTri()
        points, edges = vtk_to_tri.generate(self.network)
        self.points = points
        self.edges = edges
        
    def inside(self, x, on_boundary):
         
        if len(x) == 2:
            position = np.array((x[0], x[1], 0.0))
        else:
            position = np.array(x)
        
        am_inside = False
        for eachEdge in self.edges:
            if len(self.points[eachEdge[0]]) == 2:
                e1loc = np.array((self.points[eachEdge[0]][0], self.points[eachEdge[0]][1], 0.0))
                e2loc = np.array((self.points[eachEdge[1]][0], self.points[eachEdge[1]][1], 0.0))
            else:
                e1loc = np.array(self.points[eachEdge[0]])
                e2loc = np.array(self.points[eachEdge[1]])
            if vtk.vtkLine.DistanceToLine(position, e1loc, e2loc) <= self.tol:
                dp1 = np.linalg.norm(e1loc - position)
                dp2 = np.linalg.norm(e2loc - position)
                dpLine = np.linalg.norm(e1loc - e2loc)
                if dp1 + dp2 <= dpLine + self.tol:
                    am_inside = True
                    break
        return am_inside