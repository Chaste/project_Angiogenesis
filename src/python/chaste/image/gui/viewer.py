import logging
import wx

import chaste.gui.properties

if chaste.gui.properties._have_wx:
    import chaste.gui.panels.base
    
    class ImageViewer2dPanel(chaste.gui.panels.base.Panel):
        
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
            self.name = "ImageViewer2d"
            
        def add_controls(self):
            
            ''' 
            Add wx controls to the panel
            '''  
            
            self.select_input_file = wx.Button(parent = self, label="Load Image File")
            self.file_label = wx.StaticText(parent = self, label="Input File: ")
            self.file_text = wx.StaticText(parent = self, label = "None")
    
        def size_controls(self):
            
            ''' 
            Size the controls
            '''  
            
            centre_fmt = [0, wx.CENTER, 3]
            expand_fmt = [1, wx.EXPAND, 3]
            vbox = wx.BoxSizer(wx.VERTICAL)
            vbox.AddSpacer(10)
            vbox.Add(self.select_input_file, *centre_fmt)
            vbox.Add(self.file_label, *centre_fmt)
            vbox.Add(self.file_text, *centre_fmt)
            vbox.AddSpacer(10)
    
            self.SetSizer(vbox)
            vbox.Fit(self) 
            
        def bind_events(self):
            
            ''' 
            Bind the events
            '''  
            
            self.select_input_file.Bind(wx.EVT_BUTTON, self.on_load_file)
            
        def on_load_file(self, event = None):
            
            self.file_name = self.get_file_name()
            self.canvas = self.GetTopLevelParent().get_2d_canvas(show = True)
            self.canvas.add_tiff(self.file_name)
             
            logging.info("Loaded Image File: " + str(self.file_name))