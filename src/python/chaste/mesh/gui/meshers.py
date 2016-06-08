import logging
import os
import wx
import dolfin as df

import chaste.utility.readwrite
import chaste.geometry.boundary_markers
import chaste.mesh.meshing2d
import chaste.mesh.converters
import chaste.plot.two.glyphs
   
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
        self.surface = chaste.utility.readwrite.read_vtk_surface(self.file_name, True, True)
        
        bound_extractor = chaste.geometry.boundary_markers.BoundaryExtractor2d()
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
        self.domain = chaste.utility.readwrite.read_vtk_surface(self.file_name)
        self.lines = chaste.interfaces.converters.other.vtkpolygon_to_lines(self.domain)
    
        glyph = chaste.plot.two.glyphs.VtkLinesGlyph(self.lines, color="blue")
        self.canvas = self.GetTopLevelParent().get_2d_canvas(show = False)
        self.canvas.add_glyph(glyph, True)
        
        if not self.canvas_set:
            self.setup_canvas()
            
        conv = chaste.mesh.converters.VtkToTri()
        points, edges = conv.generate(self.lines)
        myedges = []
        for eachEdge in edges:
            myedges.append([points[eachEdge[0]], points[eachEdge[1]]])
        self.domain_edges = myedges
         
        logging.info("Loaded Domain File: " + str(self.file_name))
        
    def on_load_centreline_file(self, event = None):
        
        self.file_name = self.get_file_name()
        self.centres = chaste.utility.readwrite.read_vtk_surface(self.file_name, True, True)
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
        
        converter = chaste.mesh.meshing2d.DolfinConverter2d()
        converter.boundary_edges = self.boundary_edges
        converter.domain_edges = self.domain_edges
        converter.network = self.network

        main_mesh, r1_mesh, r2_mesh = converter.update(self.mesh, regions_set=self.tool.region1_points is not None)
        
        if converter.new_network is not None:
            file_mod = os.path.splitext(self.file_name)[0] + "_new_network.vtp"
            chaste.utility.readwrite.write(converter.new_network, file_mod)
        
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
        
        self.tool = chaste.mesh.meshing2d.Mesher2d()
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
