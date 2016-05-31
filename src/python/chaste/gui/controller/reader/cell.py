import casie.gui.properties
import casie.controller.base
import casie.rwc.vtk_reader

if casie.gui.properties._is_gui_session:
    import matplotlib
    matplotlib.use('WXAgg')
    import matplotlib.pyplot as mpp

class Controller(casie.controller.base.Controller):
    
    def __init__(self, controller_id = 0, label = None, io_data = None):
        
        casie.controller.base.Controller.__init__(self, controller_id, "reader", "cell", 
                                                  label, io_data)

        self.data = None
        self.file_path = None

    def activate(self):
        success = False
        if self.file_path is not None:
            self.data = self.read_from_file()
            success = True
            
        if casie.gui.properties._is_gui_session:
            if success:
                self.glyph.status_button.set_color("green")
            else:
                self.glyph.status_button.set_color("red")
            
        return self.own_type, self.data  
    
    def set_file_path(self, filePath):
        self.file_path = filePath  
        
    def read_from_file(self):
        reader = casie.rwc.vtk_reader.CellReader()
        return reader.read(self.file_path)
            
if casie.gui.properties._is_gui_session:
    import wx
    
    class Panel(casie.controller.base.Panel):
    
        ''' Panel for holding the configuration options and canvas panels
        '''
        
        def __init__(self, parent):
            casie.controller.base.Panel.__init__(self, parent, "reader", "cell")
    
        def add_gui_controls(self):
            self.cell_reader_button = wx.Button(parent = self, label="Choose File")
            self.cell_reader_button.Bind(wx.EVT_BUTTON, self.on_choose_file)
            self.sizer.Add(self.cell_reader_button, 0, wx.CENTER)
            
        def on_view_2d(self, event):
            if self.active_controller is not None:
                if self.active_controller.data is None:
                    self.active_controller.activate()
                fig = mpp.figure(facecolor="white") 
                fig.ax = fig.add_subplot(111) 
                try:
                    fig.glyphs = self.active_controller.data.get_glyphs()
                    for eachGlyph in fig.glyphs:
                        eachGlyph.attach_to_axes(fig.ax)
                except:
                    print "No glyphs defined for this controller"
                fig.ax.autoscale_view()
                self.GetTopLevelParent().set_scene_2d(fig)
            
        def on_choose_file(self, event):
            file_dialogue = wx.FileDialog(self, style=wx.FD_OPEN)
            file_dialogue.ShowModal()
            
            if self.active_controller is not None:
                self.active_controller.set_file_path(file_dialogue.GetPath())
            file_dialogue.Destroy()