import casie.gui.properties
import casie.controller.base

import casie.population.vessel

if casie.gui.properties._is_gui_session:
    import matplotlib
    matplotlib.use('WXAgg')
    import matplotlib.pyplot as mpp
    
class Controller(casie.controller.base.Controller):
    
    def __init__(self, controller_id = 0, label = None, io_data = None):
        
        casie.controller.base.Controller.__init__(self, controller_id, "generator", "vessel", 
                                                  label, io_data)

        self.data = None
        self.status = "deactivated"

    def activate(self):
        self.status = "error"
        network = casie.population.vessel.VesselNetwork()
        
        for eachConnector in self.input_data_connectors:
            controller_type, data = eachConnector.activate()
            if (controller_type[0] == "generator" or controller_type[0] == "filter") and controller_type[1] == "part":
                vessels = self.generate_vessels(data)
                network.add_vessels(vessels)
                self.data = network
                self.status = "activated"
                
        if casie.gui.properties._is_gui_session:
            if self.status == "activated":
                self.glyph.status_button.set_color("green")
            else:
                self.glyph.status_button.set_color("red")
                
        return self.own_type, self.data  
    
    def generate_vessels(self, part):
        
        part.set_points_ids()
        nodes = []
        for eachPoint in part.points:
            nodes.append(casie.population.vessel.VesselNode(eachPoint.location))
        
        vessels = []
        for eachEdge in part.edges:
            segment = casie.population.vessel.VesselSegment(nodes[eachEdge.start_point.id], nodes[eachEdge.end_point.id])
            vessel = casie.population.vessel.Vessel([segment,])
            vessels.append(vessel)
        return vessels 
    
if casie.gui.properties._is_gui_session:
    
    class Panel(casie.controller.base.Panel):
    
        ''' Panel for holding the configuration options and canvas panels
        '''
        
        def __init__(self, parent):
            casie.controller.base.Panel.__init__(self, parent, "generator", "vessel")
    
        def add_gui_controls(self):
            pass
        
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
                self.GetTopLevelParent().canvas_manager.set_scene_2d_canvas(fig)