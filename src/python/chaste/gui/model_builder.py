''' 2D Wx-matplotlib Graphics Window
'''

import wx

import casie.gui.panels.base

class ModelBuilderPanel(casie.gui.panels.base.Panel):
        
    def __init__(self, parent):
        
        ''' 
        Set up the panel, add the controls
        '''
        
        casie.gui.panels.base.Panel.__init__(self, parent)
        self.name = "ModelBuilder"
        
    def add_controls(self):
        
        ''' 
        Add wx controls to the panel
        '''  
        
        self.show = wx.Button(parent = self, label="Show")
        
    def size_controls(self):
        
        ''' 
        Size the controls
        '''  
        
        centre_fmt = [0, wx.CENTER, 3]
        expand_fmt = [1, wx.EXPAND, 3]
        vbox = wx.BoxSizer(wx.VERTICAL)
        vbox.AddSpacer(10)
        vbox.Add(self.show, *centre_fmt)
        vbox.AddSpacer(10)

        self.SetSizer(vbox)
        vbox.Fit(self) 
        
    def bind_events(self):
        
        ''' 
        Bind the events
        '''  
        
        self.show.Bind(wx.EVT_BUTTON, self.on_show)
        
    def on_show(self, event = None):
        
        self.canvas = self.GetTopLevelParent().get_mb_canvas(show = True)