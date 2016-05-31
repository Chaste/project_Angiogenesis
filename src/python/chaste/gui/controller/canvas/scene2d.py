import casie.gui.properties
import casie.controller.base

import matplotlib
matplotlib.use('WXAgg')
import matplotlib.pyplot as mpp

class Controller(casie.controller.base.Controller):
    
    def __init__(self, controller_id = 0, label = None, io_data = None):
        
        casie.controller.base.Controller.__init__(self, controller_id, "canvas", "scene2d", 
                                                  label, io_data)

        self.file_path = None

    def activate(self):
        # Set up a blank figure
        fig = mpp.figure(facecolor="white") 
        fig.ax = fig.add_subplot(111) 
        fig.glyphs = []

        for eachConnector in self.input_data_connectors:
            controller_type, data = eachConnector.activate()

            if controller_type[0] == "reader" or "generator": 
                if controller_type[1] == "cell":
                    glyphs = data.get_glyphs()
                    for eachGlyph in glyphs:
                        eachGlyph.attach_to_axes(fig.ax)
                    fig.glyphs.append(glyphs)
                if controller_type[1] == "mesh":
                    collection = data.get_glyphs().collection
                    fig.ax.add_collection(collection) 
                if controller_type[1] == "vessel":
                    glyphs = data.get_glyphs()
                    for eachGlyph in glyphs:
                        eachGlyph.attach_to_axes(fig.ax)
                    fig.glyphs.append(glyphs)                     
                fig.ax.autoscale_view()
        return fig       
    
if casie.gui.properties._is_gui_session:
    
    class Panel(casie.controller.base.Panel):
    
        ''' Panel for holding the configuration options and canvas panels
        '''
        
        def __init__(self, parent):
            casie.controller.base.Panel.__init__(self, parent, "canvas", "scene2d")
    
        def add_gui_controls(self):
            pass
        
        def on_view_2d(self, event):
            if self.active_controller is not None:
                fig = self.active_controller.activate()
                self.get_canvas_manager().set_scene_2d_canvas(fig)
                