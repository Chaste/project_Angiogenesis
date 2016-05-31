import casie.gui.properties
import casie.controller.base

import casie.geometry.part

if casie.gui.properties._is_gui_session:
    import matplotlib
    matplotlib.use('WXAgg')
    import matplotlib.pyplot as mpp

class Controller(casie.controller.base.Controller):
    
    def __init__(self, controller_id = 0, label = None, io_data = None):
        
        casie.controller.base.Controller.__init__(self, controller_id, "generator", "part", 
                                                  label, io_data)

        self.data = None

    def activate(self):
        return self.own_type, self.data 
    
    def add_part(self, part):
        if self.data is None:
            self.data = part
        else:
            self.data.append(part)
    
    def set_part(self, part):
        self.data = part 
    
if casie.gui.properties._is_gui_session:
    import wx
    
    class Panel(casie.controller.base.Panel):
    
        ''' Panel for holding the configuration options and canvas panels
        '''
        
        def __init__(self, parent):
            casie.controller.base.Panel.__init__(self, parent, "generator", "part")
    
        def add_gui_controls(self):
            self.start_sketch_button = wx.Button(parent = self, label="Sketch")
            self.start_sketch_button.Bind(wx.EVT_BUTTON, self.on_sketch)
            
            self.square_button = wx.Button(parent = self, label="Add Square")
            self.square_button.Bind(wx.EVT_BUTTON, self.on_square)
            
            self.circle_button = wx.Button(parent = self, label="Add Circle")
            self.circle_button.Bind(wx.EVT_BUTTON, self.on_circle)
            
            self.sizer.Add(self.start_sketch_button, 0, wx.CENTER)
            self.sizer.Add(self.square_button, 0, wx.CENTER)
            self.sizer.Add(self.circle_button, 0, wx.CENTRE)
            
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
            
        def on_sketch(self, event):
            self.GetTopLevelParent().panel.scene.show_2d()
            self.GetTopLevelParent().panel.scene.scene2d_canvas.start_temp_sketch(self)
            
        def on_square(self, event):
            if self.active_controller is not None:
                generator = casie.geometry.part.Generator()
                part = generator.make_square()
                self.active_controller.add_part(part)
                
        def on_circle(self, event):
            if self.active_controller is not None:
                generator = casie.geometry.part.Generator()
                part = generator.make_circle()
                self.active_controller.add_part(part)
            