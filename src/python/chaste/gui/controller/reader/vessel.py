import casie.gui.properties
import casie.controller.base
#import casie.reader_writer.vtk_reader


class Controller(casie.controller.base.Controller):
    
    def __init__(self, controller_id = 0, label = None, io_data = None):
        
        casie.controller.base.Controller.__init__(self, controller_id, "reader", "vessel", 
                                                  label, io_data)

        self.file_path = None

    def activate(self):
        
        output_types = []
        output_types.append(self.own_type)
        
        output_data = []
            
        if self.file_path is not None:
            own_data = self.read_from_file()
        else:
            own_data = None
        
        output_data.append(own_data)  
            
        return output_types, output_data  
    
    def set_file_path(self, filePath):
        self.file_path = filePath  
        
    def read_from_file(self):
        pass
        return 0.0
    
if casie.gui.properties._is_gui_session:
    import wx
    
    class Panel(casie.controller.base.Panel):
    
        ''' Panel for holding the configuration options and canvas panels
        '''
        
        def __init__(self, parent):
            casie.controller.base.Panel.__init__(self, parent, "reader", "vessel")
    
        def add_gui_controls(self):
            self.cell_reader_button = wx.Button(parent = self, label="Choose File")
            
            self.cell_reader_button.Bind(wx.EVT_BUTTON, self.on_choose_file)
            self.sizer.Add(self.cell_reader_button, 1, wx.EXPAND)
            
        def on_choose_file(self, event):
            file_dialogue = wx.FileDialog(self, style=wx.FD_OPEN)
            file_dialogue.ShowModal()
            
            if self.active_controller is not None:
                self.active_controller.set_file_path(file_dialogue.GetPath())
            file_dialogue.Destroy()