import casie.gui.properties
import casie.controller.base


class Controller(casie.controller.base.Controller):
    
    def __init__(self, controller_id = 0, label = None, io_data = None):
        
        casie.controller.base.Controller.__init__(self, controller_id, "simulation", "discrete", 
                                                  label, io_data)

    def activate(self):
        pass 
    
if casie.gui.properties._is_gui_session:
    
    import wx
    
    class Panel(casie.controller.base.Panel):
    
        ''' Panel for holding the configuration options and canvas panels
        '''
        
        def __init__(self, parent):
            casie.controller.base.Panel.__init__(self, parent, "simulation", "discrete")
    
        def add_gui_controls(self):
            
            pass
#             self.cell_consumption_rate_label = wx.StaticText(parent = self, label="Cell Consumption Rate")
#             self.cell_consumption_rate_ctrl = wx.TextCtrl(self, value="0.001")
#             
#             
#             
#             self.value_label = wx.StaticText(parent = self, label="Value")
#             self.prop1_name_ctrl = wx.TextCtrl(self)
#             self.prop1_value_ctrl = wx.TextCtrl(self)
#             self.prop2_name_ctrl = wx.TextCtrl(self)
#             self.prop2_value_ctrl = wx.TextCtrl(self)  
#             
#             flex_grid = wx.FlexGridSizer(3, cols=2, vgap=3, hgap=3)    
#             flex_grid.Add(self.name_label, 1, wx.EXPAND) 
#             flex_grid.Add(self.value_label, 1, wx.EXPAND)
#             flex_grid.Add(self.prop1_name_ctrl, 1, wx.EXPAND) 
#             flex_grid.Add(self.prop1_value_ctrl, 1, wx.EXPAND)
#             flex_grid.Add(self.prop2_name_ctrl, 1, wx.EXPAND) 
#             flex_grid.Add(self.prop2_value_ctrl, 1, wx.EXPAND)           
#             
#             self.set_prop_button = wx.Button(parent = self, label="Set")
#             self.set_prop_button.Bind(wx.EVT_BUTTON, self.on_set_props)
#             self.sizer.Add(flex_grid, 1, wx.EXPAND)
#             self.sizer.Add(self.set_prop_button, 1, wx.CENTER)