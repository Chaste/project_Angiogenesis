''' 2D Wx-matplotlib Graphics Window
'''
import random

import wx
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas

import casie.gui.controller.connector


class Panel(wx.Panel):
 
    ''' Matplotlib panel for the model builder
    '''   
    
    def __init__(self, parent):
        wx.Panel.__init__(self, parent)
        self.parent = parent
        self.dpi = 100
        self.color = "white"
        self.SetMaxSize((1200, 800))
        self.figure = Figure(facecolor=self.color, dpi = self.dpi)
        self.free_connector = None
        self.viewer_connection = None
        self.controller_glyphs = []
        self.connector_glyphs = []
        self.initialize_canvas()
        
    def initialize_canvas(self):
        self.width = 100
        self.height = 100
        self.figure.canvas = FigureCanvas(self, -1, self.figure)  
        self.ax = self.figure.add_subplot(111) 
        self.ax.get_xaxis().set_ticks([])
        self.ax.get_yaxis().set_ticks([])
        
        self.ax.set_xlim(0, self.width)
        self.ax.set_ylim(0, self.height)
        self.figure.subplots_adjust(left=0.02, bottom=0.02, right=0.98, top=0.98)
        
        self.figure.canvas.mpl_connect('key_press_event', self.on_key_press)
        self.figure.canvas.mpl_connect('motion_notify_event', self.on_motion)
        self.figure.canvas.mpl_connect('pick_event', self.on_pick)
        self.figure.canvas.mpl_connect('button_release_event', self.on_release)
            
        self.vbox = wx.BoxSizer(wx.VERTICAL)
        self.vbox.Add(self.figure.canvas, 1, wx.LEFT | wx.TOP | wx.EXPAND)
        self.SetSizer(self.vbox)
        self.Fit()
        self.SendSizeEvent()
        
    def add_controller_glyph(self, glyph):
        position = [0.1*self.width + 0.8*random.random()*self.width, 
                    0.1*self.height + 0.8*random.random()*self.height]
        glyph.move_to(position)
        glyph.attach_to_axes(self.ax)
        self.controller_glyphs.append(glyph)
        self.re_draw()  
        
    def add_connector_glyph(self, glyph):
        self.connector_glyphs.append(glyph)
        glyph.attach_to_axes(self.ax)
        
    def get_controller_manager(self):
        return self.GetTopLevelParent().controller_manager
    
    def remove_controller_glyph(self, glyph):
        self.controller_glyphs.remove(glyph)
        glyph.remove()
        
    def remove_connector_glyph(self, glyph):
        self.connector_glyphs.remove(glyph)
        glyph.line.remove()
        
    def re_draw(self):
        self.ax.autoscale_view()
        self.figure.canvas.draw()  
        
    def on_key_press(self, event):
        if self.free_connector is not None and "Esc" in event.key:
            self.free_connector.remove()
            self.free_connector = None
            self.re_draw()

    def on_motion(self, event):
        if self.free_connector is not None:
            self.free_connector.set_end_point([event.xdata, event.ydata])
            self.re_draw() 
        else:
            for eachControllerGlyph in self.controller_glyphs:
                if eachControllerGlyph.pressed:
                    eachControllerGlyph.move_to([event.xdata, event.ydata])
                    self.re_draw()   
        
    def on_release(self, event):
        for eachControllerGlyph in self.controller_glyphs:
            eachControllerGlyph.pressed = False
        self.re_draw()  
    
    def on_pick(self, event):
        for eachControllerGlyph in self.controller_glyphs:     
            if event.artist == eachControllerGlyph.output_port:
                if self.free_connector is None:
                    self.free_connector = casie.gui.controller.connector.Glyph()
                    self.free_connector.attach_to_axes(self.ax)
                    self.free_connector.set_start_controller_glyph(eachControllerGlyph)
                    self.re_draw() 
                    break
            elif event.artist == eachControllerGlyph.input_port and self.free_connector is not None:
                self.get_controller_manager().new_connector(self.free_connector.start_controller.controller, 
                                                                   eachControllerGlyph.controller)
                self.free_connector.line.remove()
                self.free_connector = None
                self.re_draw() 
                break          
            elif event.artist == eachControllerGlyph.rect:   
                if event.mouseevent.dblclick:
                    eachControllerGlyph.controller.panel.set_active_controller(eachControllerGlyph.controller)
                    self.GetTopLevelParent().set_config_panel(eachControllerGlyph.controller.panel)
                elif event.mouseevent.key is not None and "control" in event.mouseevent.key:
                    self.get_controller_manager().remove_contoller(eachControllerGlyph.controller)
                    self.remove_controller_glyph(eachControllerGlyph)
                else:    
                    eachControllerGlyph.pressed = True
                break
               
        for eachConnectorGlyph in self.connector_glyphs:
            if event.artist == eachConnectorGlyph.line:
                if event.mouseevent.key is not None and "control" in event.mouseevent.key:
                    self.get_controller_manager().remove_connector(eachConnectorGlyph.connector)
                    self.remove_connector_glyph(eachConnectorGlyph)