import wx

import matplotlib
matplotlib.use("wxAgg")
import casie.gui.panels.base
import casie.gui.sub_panels.canvas2d
import casie.gui.sub_panels.canvas3d
import casie.gui.sub_panels.model_builder_canvas

class Panel(casie.gui.panels.base.Panel):
    
    ''' 
    Panel for displaying 2d and 3d canvases and their options controls
    
    Attributes:
    -----------
    
    scene2d_canvas : a 2D matplotlib canvas
    scene3d_canvas : a 3D VTK render window
    '''
    
    def __init__(self, parent):
        
        ''' 
        Set up the panel, add the controls, bind the events
        '''
        
        casie.gui.panels.base.Panel.__init__(self, parent)
        
    def add_controls(self):
        
        ''' 
        Add wx controls to the panel
        '''  
        
        self.canvas2d = casie.gui.sub_panels.canvas2d.Panel(self)
        self.canvas3d = casie.gui.sub_panels.canvas3d.Panel(self)
        self.mbcanvas = casie.gui.sub_panels.model_builder_canvas.Panel(self)
        self.canvas3d.Hide()
        self.mbcanvas.Hide()
        
    def size_controls(self):
        
        ''' 
        Size the controls
        ''' 
        
        vbox = wx.BoxSizer(wx.VERTICAL)
        vbox.Add(self.canvas2d, 6, wx.EXPAND)
        vbox.Add(self.canvas3d, 6, wx.EXPAND)
        vbox.Add(self.mbcanvas, 6, wx.EXPAND)
        self.SetSizer(vbox)
        vbox.Fit(self) 
        
    def show_3d(self, event = None):
        
        ''' 
        Hide the 2D canvas and show the 3D one
        ''' 
        
        self.canvas3d.Show()
        self.canvas2d.Hide()
        self.mbcanvas.Hide()
        self.SendSizeEvent()
        
    def show_2d(self, event = None):
        
        ''' 
        Hide the 3D canvas and show the 2D one
        ''' 
        
        self.mbcanvas.Hide()
        self.canvas3d.Hide()
        self.canvas2d.Show()
        self.SendSizeEvent()
        
    def show_model_builder(self, event = None):
        
        ''' 
        Hide the 3D canvas and show the 2D one
        ''' 
        
        self.mbcanvas.Show()
        self.canvas3d.Hide()
        self.canvas2d.Hide()
        self.SendSizeEvent()