import os
import logging
import wx
import chaste.utility
import chaste.plot.two.glyphs
import chaste.gui.panels.base
import chaste.geometry.centrelines2d

class Centrelines2dPanel(chaste.gui.panels.base.Panel):
    
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
        self.surface = chaste.utility.input_output.read_vtk_surface(self.file_name)
 
        glyph = chaste.visualization.two.glyphs.VtkLinesGlyph(self.surface)
        
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
        
        self.tool = chaste.geometry.centrelines2d.Centrelines2d()
        
    def run_tool(self):
        
        self.tool.set_surface(self.surface)
        self.tool.update()
        self.network = self.tool.get_output()
        
        glyph = chaste.visualization.two.glyphs.VesselNetworkGlyph(self.network)
        self.canvas.add_glyph(glyph, clear = False)
        
        logging.info("Extracted Centreline")