''' The Casie GUI app and main frame.

'''

import os
import logging
import wx
import casie.gui.properties
import casie.gui.menus
import casie.gui.panels.logger
import casie.gui.panels.session
import casie.gui.panels.scene
import casie.gui.panels.tools

class App(wx.App):
    
    ''' The Casie GUI app. Inherits from a wx App.
    
    Attributes
    ----------
    frame : The main GUI frame (or window in non-wx speak)
    
    '''
    
    def __init__(self):
        
        ''' Set up the main GUI frame
        
        '''
        
        wx.App.__init__(self, redirect=False)
        self.frame = Frame(None, "Casie: Cancer Simulation Environment")
        
    def launch(self):
        
        ''' Start the event handling loop
        
        '''
        
        self.MainLoop()

class Frame(wx.Frame):
    
    ''' The main GUI frame
    
    Attributes
    ----------
    session_manager : Event handling for sessions
    controller_manager : Event handling for controllers and connectors
    panel : The main panel for the frame, mostly for consistent cross-platform coloring.
    panel.session : GUI controls for managing sessions
    panel.model_builder : GUI controls for building new models
    panel.scene : GUI controls for 2D and 3D plotting and sketching
    '''
    
    def __init__(self, parent, title):
        
        ''' Set up the main frame.
        
        Create managers for event handling and the main menus and panels.
        
        ''' 
        
        wx.Frame.__init__(self, parent, title=title)
        self.SetBackgroundColour(casie.gui.properties._colors["background"]) 
        self.work_dir = os.getcwd()
        self.create_menus()
        self.create_panels()
        self.Show(True)
                
    def create_menus(self):
        
        ''' 
        Create the main menu bar for the GUI and bind events.
        '''      
        
        self.menu_bar = casie.gui.menus.MainMenuBar()
        self.SetMenuBar(self.menu_bar)
        self.Bind(wx.EVT_MENU, self.exit, self.menu_bar.file_menu.exit)
    
    def create_panels(self):
        
        '''
         Create and size the main panels
        ''' 
        
        self.panel = wx.Panel(self, wx.ID_ANY)
        self.panel.session = casie.gui.panels.session.Panel(self.panel)
        self.panel.tools = casie.gui.panels.tools.Panel(self.panel)
        self.panel.scene = casie.gui.panels.scene.Panel(self.panel)
        self.panel.logging = casie.gui.panels.logger.Panel(self.panel)
        
        left_vbox = wx.BoxSizer(wx.VERTICAL)
        left_vbox.Add(self.panel.session, 1, wx.EXPAND)
        left_vbox.Add(self.panel.tools, 3, wx.EXPAND)
        
        right_vbox = wx.BoxSizer(wx.VERTICAL)
        right_vbox.Add(self.panel.scene, 5, wx.EXPAND)
        right_vbox.Add(self.panel.logging, 1, wx.EXPAND)
        
        hbox = wx.BoxSizer()
        hbox.Add(left_vbox, 1, wx.EXPAND)
        hbox.Add(right_vbox, 3, wx.EXPAND)
        self.panel.SetSizer(hbox)   
        
        frame_sizer = wx.BoxSizer()
        frame_sizer.Add(self.panel, 1, wx.EXPAND)
        frame_sizer.Fit(self)
        
        logging.info("Launched GUI with working dir:" + self.work_dir)
        
    def exit(self, event):
        
        ''' Destroy the frame. This should also end the app.
        ''' 
        
        self.Destroy()
    
    def get_2d_canvas(self, show=False):
        if show:
            self.panel.scene.show_2d()
        return self.panel.scene.canvas2d
    
    def get_3d_canvas(self, show=False):
        if show:
            self.panel.scene.show_3d()
        return self.panel.scene.canvas3d
    
    def get_mb_canvas(self, show=False):
        if show:
            self.panel.scene.show_model_builder()
        return self.panel.scene.mbcanvas
    
    def change_tool(self, name):
        new_panel = self.panel.tools.get_panel(name)
        self.panel.tools.show_panel(new_panel)
        
if __name__ == "__main__":

    
    # Launch the window
    app = App()
    app.launch()