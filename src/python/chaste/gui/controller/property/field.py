import casie.gui.properties
import casie.controller.base


class Controller(casie.controller.base.Controller):
    
    def __init__(self, controller_id = 0, label = None, io_data = None):
        
        casie.controller.base.Controller.__init__(self, controller_id, "property", "field", 
                                                  label, io_data)
        
        self.props = None
        
    def activate(self):
        return self.own_type, self.props  
    
    def set_props(self, props):
        self.props = props  
    
if casie.gui.properties._is_gui_session:
    import wx
    
    class Panel(casie.controller.base.Panel):
    
        ''' Panel for holding the configuration options and canvas panels
        '''
        
        def __init__(self, parent):
            casie.controller.base.Panel.__init__(self, parent, "property", "field")
    
        def add_gui_controls(self):
            self.consumption_rate_label = wx.StaticText(parent = self, label="Name")
            self.consumption_rate_ctrl = wx.TextCtrl(self, value="Oxygen")
            self.production_rate_label = wx.StaticText(parent = self, label="Diffusivity")
            self.production_rate_ctrl = wx.TextCtrl(self, value="1.0")            
            
            flex_grid = wx.FlexGridSizer(2, cols=2, vgap=3, hgap=3)    
            flex_grid.Add(self.consumption_rate_label, 1, wx.EXPAND) 
            flex_grid.Add(self.consumption_rate_ctrl, 1, wx.EXPAND)
            flex_grid.Add(self.production_rate_label, 1, wx.EXPAND) 
            flex_grid.Add(self.production_rate_ctrl, 1, wx.EXPAND)          
            
            self.set_prop_button = wx.Button(parent = self, label="Set")
            self.set_prop_button.Bind(wx.EVT_BUTTON, self.on_set)
            self.sizer.Add(flex_grid, 1, wx.EXPAND)
            self.sizer.Add(self.set_prop_button, 1, wx.CENTER)
            
        def on_set(self, event):

            if self.active_controller is not None:
                props = {"Name" :self.consumption_rate_ctrl.GetValue(), 
                         "Diffusivity" : float(self.production_rate_ctrl.GetValue())}
                self.active_controller.set_props(props)