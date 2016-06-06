''' 2D Wx-matplotlib Graphics Window
'''

import wx
import numpy as np
import matplotlib.pyplot as mpp
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas

import chaste.gui.sketcher.sketcher

class Panel(wx.Panel):
 
    ''' Matplotlib panel for 2D plotting
    '''   
    
    def __init__(self, parent):
        wx.Panel.__init__(self, parent)
        
        self.parent = parent
        self.dpi = 100
        self.color = "white"
        
        self.figure = Figure(facecolor=self.color, dpi = self.dpi)
        self.canvas = FigureCanvas(self, -1, self.figure)
        
        self.clear()
        
        self.add_grid()
        
        self.start_temp_sketch(self)
        
        self.vbox = wx.BoxSizer(wx.VERTICAL)
        self.vbox.Add(self.canvas, 1, wx.LEFT | wx.TOP | wx.EXPAND)
        self.SetSizer(self.vbox)
        self.Fit() 
        
    def clear(self):
        
        self.figure.clear()
        self.figure.ax = self.figure.add_subplot(111) 
        
        self.figure.ax.set_xlim(0, 100)
        self.figure.ax.set_ylim(0, 100)
        self.figure.ax.set_aspect('equal', adjustable='box')
        
        self.canvas.draw()
        self.SendSizeEvent()
        
    def add_grid(self):
        
        grid_spacing = 10.0
        xlims = self.figure.ax.get_xlim()
        xlen = xlims[1] - xlims[0]
        ylims = self.figure.ax.get_ylim()
        ylen = ylims[1] - ylims[0]
        
        Xtick = np.linspace(xlims[0], xlims[1], int(xlen / grid_spacing))
        self.figure.ax.set_xticks(Xtick)
        for eachPoint in Xtick:
            self.figure.ax.axvline(x=eachPoint,ls='-', color='black')
          
        Ytick = np.linspace(ylims[0], ylims[1], int(ylen / grid_spacing))
        self.figure.ax.set_yticks(Ytick)
        for eachPoint in Ytick:
            self.figure.ax.axhline(y=eachPoint,ls='-', color='black')
        
    def add_image(self, figure):
        
        self.figure = figure
        self.canvas.figure = self.figure
        self.canvas.draw()
        self.SendSizeEvent()
        
    def add_tiff(self, path):
        
        self.clear()
        
        self.figure.ax.imshow(mpp.imread(path))
        self.figure.ax.autoscale()
        
#        self.add_grid()
        self.start_temp_sketch(self)
            
        self.canvas.draw()
        self.SendSizeEvent()
        
    def add_glyph(self, glyph, clear = False):
        
        if clear:
            self.clear()
            
        glyph.attach_to_axes(self.figure.ax)
        self.figure.ax.autoscale()
        
        xlims = self.figure.ax.get_xlim()
        xlen = xlims[1] - xlims[0]
        ylims = self.figure.ax.get_ylim()
        ylen = ylims[1] - ylims[0]
        
        self.figure.ax.set_xlim(xlims[0]-0.05*xlen, xlims[1]+0.05*xlen)
        self.figure.ax.set_ylim(ylims[0]-0.05*ylen, ylims[1]+0.05*ylen)
        self.figure.ax.set_aspect('equal', adjustable='box')
        
        if clear:
#            self.add_grid()
            self.start_temp_sketch(self)
        
        self.canvas.draw()
        self.SendSizeEvent()
        
    def re_draw(self):
        self.figure.ax.autoscale()
        self.canvas.draw()
        
    def start_temp_sketch(self, panel):
        
        # Set up the temp sketch
        self.temp_sketch = chaste.gui.sketcher.sketcher.Sketch(self.figure.ax, sketch_type = "Temp")
        self.temp_sketch.panel = panel
        self.re_draw()