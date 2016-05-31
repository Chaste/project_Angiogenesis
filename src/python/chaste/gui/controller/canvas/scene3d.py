import casie.gui.properties
import casie.controller.base

import matplotlib
matplotlib.use('WXAgg')
import matplotlib.pyplot as mpp

class Controller(casie.controller.base.Controller):
    
    def __init__(self, controller_id = 0, label = None, io_data = None):
        
        casie.controller.base.Controller.__init__(self, controller_id, "canvas", "scene3d", 
                                                  label, io_data)

        self.file_path = None

    def activate(self):
        
        # Set up a blank figure
        fig = self.set_up_scene()

        for eachConnector in self.input_data_connectors:
            controller_type, data = eachConnector.activate()

            if controller_type[0] == "reader" and controller_type[1] == "cell":
                x_locations = []
                y_locations = []
                colors = []
                for eachCell in data.cells:
                    x_locations.append(eachCell.location[0])
                    y_locations.append(eachCell.location[1])
                    colors.append(eachCell.color)
                fig.ax1.scatter(x_locations, y_locations, marker='o', color=colors) 
        
        return fig      
        
    def set_up_scene(self):

        fig = mpp.figure(facecolor="white") 
        fig.ax1 = fig.add_subplot(111) 

        return fig  
    
if casie.gui.properties._is_gui_session:
    
    class Panel(casie.controller.base.Panel):
    
        ''' Panel for holding the configuration options and canvas panels
        '''
        
        def __init__(self, parent):
            casie.controller.base.Panel.__init__(self, parent, "canvas", "scene3d")
    
        def add_gui_controls(self):
            pass