import logging
import os
import wx

import chaste.utility.readwrite
import chaste.interfaces.vtk_tools.glyphs
import chaste.gui.panels.base
import chaste.image.image_to_surface2d

class VtkImageToPolyData2dPanel(chaste.gui.panels.base.Panel):
    
    ''' 
    Default panel
    '''  
    
    def __init__(self, parent):
        
        ''' 
        Set up the panel, add the controls
        '''
        
        chaste.gui.panels.base.Panel.__init__(self, parent)
        self.name = "VtkImageToPolyData2d"
        
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
        chaste.utility.readwrite.write(self.surface, file_mod)
        
        file_mod = os.path.splitext(self.file_name)[0] + "_boundaries.vtp"
        chaste.utility.readwrite.write(self.boundaries, file_mod)
        
        logging.info("Saved Surface File to : " + str(file_mod))
        
    def setup_tool(self):
        
        self.tool = chaste.image.image_to_surface2d.VtkImageToPolyData2d()
        self.tool.input = chaste.utility.readwrite.read(self.file_name)
        
    def run_tool(self):
        
        self.tool.update()
        self.surface = self.tool.surface
        self.boundaries = self.tool.output
        
        glyph = chaste.interfaces.vtk_tools.glyphs.VtkLinesGlyph(self.boundaries)
        self.canvas.add_glyph(glyph, clear = False)
        
        logging.info("Extracted Surfaces")