import casie.gui.properties
import casie.controller.base

import casie.rwc.vtk_writer

class Controller(casie.controller.base.Controller):
    
    def __init__(self, controller_id = 0, label = None, io_data = None):
        
        casie.controller.base.Controller.__init__(self, controller_id, "writer", "mesh", 
                                                  label, io_data)

        self.file_directory = None

    def activate(self):
        
        for eachConnector in self.input_data_connectors:
            controller_type, data = eachConnector.activate()
            if controller_type[0] == "generator" and controller_type[1] == "mesh":
                if self.file_directory is not None:
                    self.write_to_file(data)
            
    def set_file_directory(self, fileDirectory):
        self.file_directory = fileDirectory  
        
    def write_to_file(self, data):
        writer = casie.rwc.vtk_writer.MeshWriter()
        writer.write(data, self.file_directory+"/mesh", 2)
    
if casie.gui.properties._is_gui_session:
    import wx
    
    class Panel(casie.controller.base.Panel):
    
        ''' Panel for holding the configuration options and canvas panels
        '''
        
        def __init__(self, parent):
            casie.controller.base.Panel.__init__(self, parent, "writer", "mesh")
    
        def add_gui_controls(self):
            self.choose_directory_button = wx.Button(parent = self, label="Choose Directory")
            self.choose_directory_button.Bind(wx.EVT_BUTTON, self.on_choose_directory)
            
            self.write_button = wx.Button(parent = self, label="Write")
            self.write_button.Bind(wx.EVT_BUTTON, self.on_write)
            
            self.sizer.Add(self.choose_directory_button, 1, wx.EXPAND)
            self.sizer.Add(self.write_button, 1, wx.EXPAND)

        def on_choose_directory(self, event):
            dir_dialogue = wx.DirDialog(self, style=wx.FD_OPEN)
            dir_dialogue.ShowModal()
            
            if self.active_controller is not None:
                self.active_controller.set_file_directory(dir_dialogue.GetPath())
            dir_dialogue.Destroy()
            
        def on_write(self, event):
            if self.active_controller is not None:
                self.active_controller.activate()          