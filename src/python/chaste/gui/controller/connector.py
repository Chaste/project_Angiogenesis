import chaste.gui
import matplotlib.lines

class Connector():
    def __init__(self, input_controller = None, output_controller = None):
        
        self.input_controller = input_controller
        self.output_controller = output_controller
        
        if self.input_controller is not None:
            self.input_controller.add_output_data_connector(self)
            
        if self.output_controller is not None:
            self.output_controller.add_input_data_connector(self)
        
        self.glyph = None
        
    def deactivate(self):
        self.input_controller.deactivate()
    
    def activate(self):
        return self.input_controller.activate()
    
    
class Glyph:

    def __init__(self, start_point = [0.0, 0.0], end_point = [1.0, 1.0]):  
        self.start_controller = None
        self.end_controller = None

        self.start_point = start_point
        self.end_point = end_point
            
        self.color = "grey"
        self.line = matplotlib.lines.Line2D([self.start_point[0], self.end_point[0]], 
                                             [self.start_point[1], self.end_point[1]], 
                                             picker=True, 
                                             color = self.color)
        
    def attach_to_axes(self, axes):
        axes.add_line(self.line)  
        
    def update_position(self):    
        new_position = [[self.start_point[0], self.end_point[0]], 
                        [self.start_point[1], self.end_point[1]]]
        self.line.set_data(new_position)
        
    def set_start_controller_glyph(self, controller_glyph):
        controller_glyph.output_connectors.append(self)
        self.start_controller = controller_glyph
        self.start_point = self.start_controller.op_position
        self.update_position()
        
    def set_end_controller_glyph(self, controller_glyph):
        controller_glyph.input_connectors.append(self)
        self.end_controller = controller_glyph
        self.end_point = self.end_controller.ip_position
        self.update_position()
        
    def set_end_point(self, position):
        self.end_point = position
        self.update_position()
        
    def set_start_point(self, position):
        self.start_point = position  
        self.update_position()