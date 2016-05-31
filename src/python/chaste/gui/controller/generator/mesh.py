import casie.gui.properties
import casie.controller.base

import casie.geometry.mesh

class Controller(casie.controller.base.Controller):
    
    def __init__(self, controller_id = 0, label = None, io_data = None):
        
        casie.controller.base.Controller.__init__(self, controller_id, "generator", "mesh", 
                                                  label, io_data)

        self.mesh = None
        
    def activate(self):
        for eachConnector in self.input_data_connectors:
            controller_type, data = eachConnector.activate()
            if controller_type[0] == "generator" and controller_type[1] == "part":
                self.generate_mesh(data)
        
        return self.own_type, self.mesh
    
    def generate_mesh(self, part):
        self.mesh = casie.geometry.mesh.Mesh()
        self.mesh.generate(part, 2)
        
    def set_mesh(self, mesh):
        self.mesh = mesh
    
if casie.gui.properties._is_gui_session:
    
    import wx

    class Panel(casie.controller.base.Panel):
    
        ''' Panel for holding the configuration options and canvas panels
        '''
        
        def __init__(self, parent):
            casie.controller.base.Panel.__init__(self, parent, "generator", "mesh")
    
        def add_gui_controls(self):
            self.reg_grid_label = wx.StaticText(parent = self, label="Regular Grid")
            self.reg_grid_extents_label = wx.StaticText(parent = self, label="Extents:")
            self.reg_grid_extents_x_label = wx.StaticText(parent = self, label="X")
            self.reg_grid_extents_x_ctrl = wx.TextCtrl(self, value="10")
            self.reg_grid_extents_y_label = wx.StaticText(parent = self, label="Y")
            self.reg_grid_extents_y_ctrl = wx.TextCtrl(self, value="10")
            self.reg_grid_extents_z_label = wx.StaticText(parent = self, label="Z")
            self.reg_grid_extents_z_ctrl = wx.TextCtrl(self, value="1")
            self.reg_grid_spacing_label = wx.StaticText(parent = self, label="Spacing:")
            self.reg_grid_spacing_ctrl = wx.TextCtrl(self, value="1")
            self.reg_grid_set_button = wx.Button(parent = self, label="Set")        
            
            flex_grid = wx.FlexGridSizer(4, cols=7, vgap=3, hgap=3)    
            flex_grid.Add(self.reg_grid_label, 1, wx.EXPAND) 
            flex_grid.Add((0,0), 1, wx.EXPAND) 
            flex_grid.Add((0,0), 1, wx.EXPAND) 
            flex_grid.Add((0,0), 1, wx.EXPAND) 
            flex_grid.Add((0,0), 1, wx.EXPAND) 
            flex_grid.Add((0,0), 1, wx.EXPAND) 
            flex_grid.Add((0,0), 1, wx.EXPAND) 
            flex_grid.Add(self.reg_grid_extents_label, 1, wx.EXPAND)
            flex_grid.Add(self.reg_grid_extents_x_label, 1, wx.EXPAND) 
            flex_grid.Add(self.reg_grid_extents_x_ctrl, 1, wx.EXPAND) 
            flex_grid.Add(self.reg_grid_extents_y_label, 1, wx.EXPAND)
            flex_grid.Add(self.reg_grid_extents_y_ctrl, 1, wx.EXPAND) 
            flex_grid.Add(self.reg_grid_extents_z_label, 1, wx.EXPAND)  
            flex_grid.Add(self.reg_grid_extents_z_ctrl, 1, wx.EXPAND)
            flex_grid.Add(self.reg_grid_spacing_label, 1, wx.EXPAND) 
            flex_grid.Add(self.reg_grid_spacing_ctrl, 1, wx.EXPAND) 
            flex_grid.Add((0,0), 1, wx.EXPAND) 
            flex_grid.Add((0,0), 1, wx.EXPAND) 
            flex_grid.Add((0,0), 1, wx.EXPAND) 
            flex_grid.Add((0,0), 1, wx.EXPAND) 
            flex_grid.Add((0,0), 1, wx.EXPAND)           
            flex_grid.Add(self.reg_grid_set_button, 1, wx.EXPAND)               
            
            self.sizer.Add(flex_grid, 1, wx.EXPAND)
            self.reg_grid_set_button.Bind(wx.EVT_BUTTON, self.on_set_regular)
            
        def on_set_regular(self, event):
            if self.active_controller is not None:
                extents = [int(self.reg_grid_extents_x_ctrl.GetValue()),
                           int(self.reg_grid_extents_y_ctrl.GetValue()),
                           int(self.reg_grid_extents_z_ctrl.GetValue()),]
                spacing = float(self.reg_grid_spacing_ctrl.GetValue())
                
                mesh = casie.geometry.mesh.RegularGrid(spacing, extents)
                self.active_controller.set_mesh(mesh)