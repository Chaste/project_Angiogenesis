import casie.gui.properties
import connector

class Controller():
    
    """ Base controller class for GUI and shell sessions
    """
    
    def __init__(self, controller_id = 0, ctr_type = "none", tartge_type = "none", 
                 name = "Controller", io_data = None):
        
        self.id = controller_id
        self.ctr_type = ctr_type
        self.target_type = tartge_type
        self.own_type = (self.ctr_type, self.target_type)
        self.name = name
        self.glyph = None
        self.panel = None
        
        self.input_data_connectors = []
        self.output_data_connectors = []
        
        if io_data is not None:
            self.generate_from_file(io_data)
    
    def add_output_data_connector(self, connector):
        self.output_data_connectors.append(connector)
        
    def add_input_data_connector(self, connector):
        self.input_data_connectors.append(connector)
        
    def add_glyph(self, glyph):
        self.glyph = glyph
        
    def add_panel(self, panel):
        self.panel = panel

    def activate(self):
        output_data = []
        for eachInputDataConnector in self.input_data_connectors:
            connector_data = eachInputDataConnector.activate()
            output_data.append(connector_data[1])
            
        return self.own_type, output_data
        
    def deactivate(self):
        
        for eachInputDataConnector in self.input_data_connectors:
            eachInputDataConnector.deactivate()
    
    def generate_from_file(self):
        pass
    
    def generate_io_data(self):
        pass

if casie.gui.properties._is_gui_session:
    import wx
    import matplotlib.patches
    import matplotlib.text

    class Panel(wx.Panel):
        
        """ Base controller panel class for the configuration options GUI panel
        """
        
        def __init__(self, parent, ctr_type = "none", target_type = "none"):
            
            wx.Panel.__init__(self, parent)
            
            self.parent = parent
            self.SetBackgroundColour(casie.gui.properties._colors[ctr_type])  
            self.SetMinSize(casie.gui.properties._sizes["ConfigOptionsCtrPanelMin"])
            
            self.ctr_type = ctr_type
            self.target_type = target_type
            self.active_controller = None
            
            self.add_controls()
            self.size_controls()
            self.bind_events()
            
        def add_controls(self):
            
            self.panel_label = wx.StaticText(parent = self, label="Target Type: " + self.target_type + 
                                             " | Controller Type: " + self.ctr_type)
            
            self.active_button = wx.Button(parent = self, label="Activate")
            self.deactivate_button = wx.Button(parent = self, label="Deactivate")
            
            self.show2d_button = wx.Button(parent = self, label="View 2D")
            self.show3d_button = wx.Button(parent = self, label="View 3D")
            
            self.controller_name_label = wx.StaticText(parent = self, label = "Name:")
            self.controller_name_ctrl = wx.TextCtrl(parent = self, value="none", style=wx.TE_PROCESS_ENTER)
            
        def size_controls(self):
            
            centre_fmt = [0, wx.CENTER, 5]
            self.label_sizer = wx.BoxSizer()
            self.label_sizer.AddSpacer(10)
            self.label_sizer.Add(self.panel_label, *centre_fmt)
            
            self.name_sizer = wx.BoxSizer()
            self.name_sizer.AddSpacer(10)
            self.name_sizer.Add(self.controller_name_label, *centre_fmt)
            self.name_sizer.Add(self.controller_name_ctrl, 1, wx.EXPAND, 5)
            self.name_sizer.Add(self.active_button, *centre_fmt)
            self.name_sizer.Add(self.deactivate_button, *centre_fmt)
            self.name_sizer.Add(self.show2d_button, *centre_fmt)
            self.name_sizer.Add(self.show3d_button,*centre_fmt)
            
            self.sizer = wx.BoxSizer(wx.VERTICAL)
            self.sizer.AddSpacer(10)
            self.sizer.Add(self.label_sizer, *centre_fmt)
            self.sizer.Add(self.name_sizer, 0, wx.LEFT)
            self.sizer.AddSpacer(10)
            self.add_gui_controls()
            self.sizer.AddSpacer(10)
    
            self.SetSizer(self.sizer)
            self.sizer.Fit(self) 
            
        def bind_events(self):

            self.active_button.Bind(wx.EVT_BUTTON, self.on_activate)
            self.deactivate_button.Bind(wx.EVT_BUTTON, self.on_deactivate)
            self.show2d_button.Bind(wx.EVT_BUTTON, self.on_view_2d)
            self.show3d_button.Bind(wx.EVT_BUTTON, self.on_view_3d)
            self.controller_name_ctrl.Bind(wx.EVT_TEXT_ENTER, self.on_controller_name_change)
            
        def get_controller_manager(self):
            
            ''' Return the controller manager
            '''       
            
            return self.GetTopLevelParent().controller_manager
        
        def get_session_manager(self):
            
            ''' Return the session manager
            '''       
            
            return self.GetTopLevelParent().session_manager
        
        def get_canvas_manager(self):
            
            ''' Return the canvas manager
            '''       
            
            return self.GetTopLevelParent().canvas_manager
            
        def on_controller_name_change(self, event):
            if self.active_controller is not None:
                self.active_controller.name = self.controller_name_ctrl.GetValue()
                if self.active_controller.glyph is not None:
                    self.active_controller.glyph.text.set_text(self.controller_name_ctrl.GetValue())
                    self.get_canvas_manager().refresh_mb_canvas()
            
        def add_gui_controls(self):
            self.label = wx.StaticText(parent = self, label=" | No Controller Selected | ")
            self.sizer.Add(self.label, 0, wx.CENTER)
            
        def set_active_controller(self, controller):
            if controller.ctr_type == self.ctr_type and controller.target_type == self.target_type:
                self.active_controller = controller
                self.controller_name_ctrl.SetValue(self.active_controller.name)
                self.update_controls()
                
        def update_controls(self):
            pass
        
        def on_view_2d(self, event):
            pass
        
        def on_view_3d(self, event):
            pass
        
        def on_activate(self, event):
            if self.active_controller is not None:
                self.active_controller.activate()
                self.get_canvas_manager().refresh_mb_canvas()
                
        def on_deactivate(self, event):
            if self.active_controller is not None:
                self.active_controller.deactivate()
                self.get_canvas_manager().refresh_mb_canvas() 
                      
    class Glyph:
        
        def __init__(self, controller, position = [0.0, 0.0]):
            
            self.controller = controller
            self.centre = position
            self.width = casie.gui.properties._sizes["MBControllerGlyphWidth"]
            self.height = casie.gui.properties._sizes["MBControllerGlyphHeight"]
            self.port_width = casie.gui.properties._sizes["MBControllerGlyphPortWidth"]
            self.port_height = casie.gui.properties._sizes["MBControllerGlyphPortHeight"]
            self.input_connectors = []
            self.output_connectors = []
            self.pressed = False
            
            self.color = [x/255.0 for x in casie.gui.properties._colors[self.controller.ctr_type][0:3]]
            self.alpha = casie.gui.properties._colors[self.controller.ctr_type][3]
            
            self.update_positions()
            self.generate_sub_glyphs()
            
        def update_positions(self): 
            self.bottom_left = [self.centre[0] - self.width/2.0, 
                                self.centre[1] - self.height/2.0]
            self.ip_position = [self.bottom_left[0], 
                                self.bottom_left[1] - self.port_height]
            self.op_position = [self.bottom_left[0] + self.width - self.port_width, 
                                self.bottom_left[1 ] - self.port_height]
            self.text_position = [self.bottom_left[0] + self.width/20.0, self.centre[1]]
            self.status_button_position = [self.centre[0] - self.port_width/2.0, 
                                          self.bottom_left[1] - self.port_height]         

        def generate_sub_glyphs(self):
            self.rect = matplotlib.patches.Rectangle(self.bottom_left, 
                                                     width = self.width, 
                                                     height = self.height, 
                                                     picker=True, 
                                                     color = self.color, 
                                                     alpha = self.alpha)
            
            self.input_port = matplotlib.patches.Rectangle(self.ip_position, 
                                                    width = self.port_width, 
                                                     height = self.port_height, 
                                                     picker=True, 
                                                     color = 'black') 
            
            self.output_port = matplotlib.patches.Rectangle(self.op_position, 
                                                    width = self.port_width, 
                                                     height = self.port_height, 
                                                     picker=True, 
                                                     color = 'black')  
            
            self.status_button = matplotlib.patches.Rectangle(self.status_button_position, 
                                                             width = self.port_width, 
                                                             height = self.port_height,
                                                             color = 'grey')    
            
            self.text = matplotlib.text.Text(x = self.text_position[0], y = self.text_position[1], 
                                             text = self.controller.name)
            
        def attach_to_axes(self, axes):
            axes.add_patch(self.rect)
            axes.add_patch(self.input_port)
            axes.add_patch(self.output_port)
            axes.add_patch(self.status_button)
            axes.add_artist(self.text)
            
        def move_sub_glyphs(self):
            self.rect.set_x(self.bottom_left[0])
            self.rect.set_y(self.bottom_left[1])
            self.input_port.set_x(self.ip_position[0])
            self.input_port.set_y(self.ip_position[1])
            self.output_port.set_x(self.op_position[0])
            self.output_port.set_y(self.op_position[1])  
            self.status_button.set_x(self.status_button_position[0])
            self.status_button.set_y(self.status_button_position[1])        
            self.text.set_position(self.text_position)
            
            for eachConnector in self.input_connectors:
                eachConnector.set_end_point(self.ip_position)
    
            for eachConnector in self.output_connectors:
                eachConnector.set_start_point(self.op_position)          
            
        def move_to(self, position):
            self.centre = position
            self.update_positions()
            self.move_sub_glyphs()
            
        def add_output_connector(self, axes):
            new_connector_glyph = connector.Glyph(axes)
            new_connector = connector.Connector()
            new_connector.input_controller = self.controller
            
            self.session.add_connector(new_connector)
            new_connector.glyph = new_connector_glyph
            new_connector_glyph.connector = new_connector
            
            new_connector_glyph.set_start_controller(self)
            self.output_connectors.append(new_connector_glyph)
            return new_connector_glyph
     
        def add_input_connector(self, connector):
            connector.set_end_controller(self)
            self.input_connectors.append(connector) 
            self.controller.add_input_data_connector(connector.connector)
            
        def on_activate(self):
            self.controller.activate()
            
        def on_deactivate(self):
            self.controller.deactivate()
            
        def remove(self):
            self.rect.remove()  
            self.input_port.remove()  
            self.output_port.remove()  
            self.status_button.remove()     
            self.text.remove()