import casie.gui.properties
import casie.controller.base

if casie.gui.properties._is_gui_session:
    import matplotlib
    matplotlib.use('WXAgg')
    import matplotlib.pyplot as mpp

class Controller(casie.controller.base.Controller):
    
    def __init__(self, controller_id = 0, label = None, io_data = None):
        
        casie.controller.base.Controller.__init__(self, controller_id, "filter", "part", 
                                                  label, io_data)

        self.data = None

    def activate(self):
        
        return self.own_type, self.data 
    
    def tesselate_input_part(self):
        
        part = None
        for eachConnector in self.input_data_connectors:
            controller_type, data = eachConnector.activate()
            if controller_type[0] == "generator" and controller_type[1] == "part":
                part = data
                break
              
        if part is not None:
            part.tesselate(merge=False)  
            self.data = part
    
if casie.gui.properties._is_gui_session:
    import wx
    
    class Panel(casie.controller.base.Panel):
    
        ''' Panel for holding the configuration options and canvas panels
        '''
        
        def __init__(self, parent):
            casie.controller.base.Panel.__init__(self, parent, "filter", "part")
    
        def add_gui_controls(self):
            self.tessellate_button = wx.Button(parent = self, label="Tessellate")
            self.tessellate_button.Bind(wx.EVT_BUTTON, self.on_tesselate)
            
            self.sizer.Add(self.tessellate_button, 0, wx.CENTER)
            
        def on_view_2d(self, event):
            if self.active_controller is not None:
                if self.active_controller.data is None:
                    self.active_controller.activate()
                fig = mpp.figure(facecolor="white") 
                fig.ax = fig.add_subplot(111) 
                if self.active_controller.data is not None:
                    fig.glyphs = self.active_controller.data.get_glyphs()
                    for eachGlyph in fig.glyphs:
                        eachGlyph.attach_to_axes(fig.ax)
                fig.ax.autoscale_view()
                self.get_canvas_manager().set_scene_2d_canvas(fig)
            
        def on_tesselate(self, event):
            if self.active_controller is not None:
                self.active_controller.tesselate_input_part()
            