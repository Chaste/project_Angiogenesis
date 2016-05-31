import casie.gui.properties
import casie.controller.base


class Controller(casie.controller.base.Controller):
    
    def __init__(self, controller_id = 0, label = None, io_data = None):
        
        casie.controller.base.Controller.__init__(self, controller_id, "writer", "cell", 
                                                  label, io_data)

        self.file_path = None

    def activate(self):
        
        if self.set_file_directory is not None and self.input_data_connectors:
            self.write_to_file(self.input_data_connectors[0].activate())
            
        return 0
    
    def set_file_directory(self, filePath):
        self.set_file_directory = filePath  
        
    def write_to_file(self, data):
        pass
    
if casie.gui.properties._is_gui_session:
    import wx
    
    class Panel(casie.controller.base.Panel):
    
        ''' Panel for holding the configuration options and canvas panels
        '''
        
        def __init__(self, parent):
            casie.controller.base.Panel.__init__(self, parent, "writer", "cell")
    
        def add_gui_controls(self):
            self.cell_reader_button = wx.Button(parent = self, label="Choose Directory")
            
            self.cell_reader_button.Bind(wx.EVT_BUTTON, self.on_choose_directory)
            self.sizer.Add(self.cell_reader_button, 1, wx.EXPAND)
            
        def on_choose_directory(self, event):
            dir_dialogue = wx.DirDialog(self, style=wx.FD_OPEN)
            dir_dialogue.ShowModal()
            
            if self.active_controller is not None:
                self.active_controller.set_file_directory(dir_dialogue.GetPath())
            dir_dialogue.Destroy()