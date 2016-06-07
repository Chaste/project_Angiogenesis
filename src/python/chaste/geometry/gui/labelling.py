import os 
import logging
import wx
import chaste.geometry.labelling
import chaste.plot.two.glyphs
import chaste.utility.rwc
import chaste.gui.panels.base

class BoundaryMarker2dPanel(chaste.gui.panels.base.Panel):
    
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
        self.name = "BoundaryMarker2d"
        self.inlets = []
        self.outlets = []
        self.select_inlets = True
        
    def add_controls(self):
        
        ''' 
        Add wx controls to the panel
        '''  
        
        self.select_input_file = wx.Button(parent = self, label="Load Boundary File")
        self.select_inlets = wx.Button(parent = self, label="Select Inlets")
        self.select_outlets = wx.Button(parent = self, label="Select Outlets")
        self.generate = wx.Button(parent = self, label="Generate")
        self.save_input_file = wx.Button(parent = self, label="Save File")

    def size_controls(self):
        
        ''' 
        Size the controls
        '''  
        
        centre_fmt = [0, wx.CENTER, 3]
        vbox = wx.BoxSizer(wx.VERTICAL)
        vbox.AddSpacer(10)
        vbox.Add(self.select_input_file, *centre_fmt)
        vbox.Add(self.select_inlets, *centre_fmt)
        vbox.Add(self.select_outlets, *centre_fmt)
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
        self.select_inlets.Bind(wx.EVT_BUTTON, self.on_select_inlets)
        self.select_outlets.Bind(wx.EVT_BUTTON, self.on_select_outlets)
        self.generate.Bind(wx.EVT_BUTTON, self.on_generate)
        self.save_input_file.Bind(wx.EVT_BUTTON, self.on_save_file)
        
    def on_select_inlets(self, event = None):
        self.select_inlets = True
        
    def on_select_outlets(self, event = None):
        self.select_inlets = False
                    
    def on_load_file(self, event = None):
        
        self.file_name = self.get_file_name()
        self.surface = chaste.utility.rwc.read_vtk_surface(self.file_name, True, triangulate=True)
 
        glyph = chaste.visualization.two.glyphs.VtkLinesGlyph(self.surface)
        
        self.canvas = self.GetTopLevelParent().get_2d_canvas(show = False)
        self.canvas.add_glyph(glyph, True)
        self.setup_canvas()
         
        logging.info("Loaded Boundary File: " + str(self.file_name))
        
    def on_save_file(self, event = None):
        
        file_mod_1 = os.path.splitext(self.file_name)[0] + "_labelled.vtp"
        file_mod_2 = os.path.splitext(self.file_name)[0] + "_open.vtp"
        chaste.utility.rwc.write_vtk_surface(file_mod_1, self.tool.get_output())
        chaste.utility.rwc.write_vtk_surface(file_mod_2, self.tool.open_surface)
        logging.info("Saved Boundary File to : " + str(file_mod_1))
        logging.info("Saved Boundary File to : " + str(file_mod_2))

    def on_generate(self, event = None):
        
        self.setup_tool()
        self.run_tool()
        
    def setup_tool(self):
        
        self.tool = chaste.geometry.labelling.BoundaryMarker2d()
        
    def run_tool(self):
        
        self.tool.set_surface(self.surface)
        self.tool.inlet_points = self.inlets
        self.tool.outlet_points = self.outlets
        self.tool.update()
        logging.info("Boundary is Labelled")

    def setup_canvas(self):
        
        self.canvas.canvas.mpl_connect('button_press_event', self.on_press)
        
    def on_press(self, event):
        
        self.canvas.temp_sketch.new_point([event.xdata, event.ydata])
        if self.select_inlets:
            self.inlets.append([event.xdata, event.ydata])
        else:
            self.outlets.append([event.xdata, event.ydata])
        self.canvas.re_draw()   