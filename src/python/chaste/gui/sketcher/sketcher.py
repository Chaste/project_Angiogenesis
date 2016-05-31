import wx
import matplotlib.lines
import matplotlib.patches

import casie.gui.panels.base

class Sketch():
    
    def __init__(self, ax, sketch_type = 'Base'):
        
        self.ax = ax
        self.type = sketch_type
        if self.type == 'Base':
            self.point_color = 'red'
            self.line_color = 'grey'
        else:
            self.point_color = 'red'
            self.line_color = 'black'          
        self.grid = None
        self.points = []
        self.lines = []
        self.arcs = []
        self.circles = []
        self.rectangles = []
        self.draw_point = True
        self.draw_line = False
        self.draw_rectangle = False
        self.draw_circle = False
        self.draw_arc = False
        
    def set_all_false(self):
        self.draw_point = False
        self.draw_line = False
        self.draw_rectangle = False
        self.draw_circle = False
        self.draw_arc = False
        
    def new_line(self, start_position, end_position = [10.0, 10.0]):
        new_line = matplotlib.lines.Line2D([start_position[0], end_position[0]], 
                                         [start_position[1], end_position[1]], 
                                         picker=True, 
                                         color = self.line_color)
        self.lines.append(new_line)
        self.ax.add_line(new_line)
        
    def set_last_line_end_point(self, position):
        last_line = self.lines[-1]
        line_data = last_line.get_data()
        line_data[0][1] = position[0]
        line_data[1][1] = position[1]
        last_line.set_data(line_data)
        
    def new_point(self, position):
        new_circle = matplotlib.patches.Circle(position,
                                               radius = 2, 
                                               color=self.point_color, 
                                               alpha = 0.5, 
                                               picker=True)
        self.points.append(new_circle)
        self.ax.add_patch(new_circle)
        
    def clear(self):
        
        for eachPoint in self.points:
            eachPoint.remove()
            
        for eachLine in self.lines:
            eachLine.remove()
            
        self.points = []
        self.lines = []
        
    def undo(self):
        
        if self.draw_point:
            self.points[-1].remove()
            del self.points[-1]
            
        if self.draw_line:
            self.lines[-1].remove()
            del self.lines[-1]
            
    def convert_to_part(self):
        
        part = casie.geometry.part.Part()
        
        points = []
        for eachPoint in self.points:
            points.append(casie.geometry.part.Point(eachPoint.center))
        part.points.extend(points)
        
        edges = []
        for eachLine in self.lines:
            for idx, eachPoint in enumerate(points):
                line_data = eachLine.get_data()
                if line_data[0][0] == eachPoint.location[0] and line_data[1][0] == eachPoint.location[1]:
                    point1 = idx
                elif line_data[0][1] == eachPoint.location[0] and line_data[1][1] == eachPoint.location[1]:
                    point2 = idx
            edges.append(casie.geometry.part.Edge(points[point1], points[point2]))
        part.edges.extend(edges)   
        
        return part
    
    def finish(self):
        
        points, edges = self.convert_to_geometry()
        self.clear()
        return points, edges
    
    def hide(self):
        for eachPoint in self.points:
            eachPoint.remove()
            
        for eachLine in self.lines:
            eachLine.remove()
            
    def show(self):
        for eachPoint in self.points:
            self.ax.add_patch(eachPoint)
            
        for eachLine in self.lines:
            self.ax.add_line(eachLine)      
    
class Panel(casie.gui.panels.base.Panel):
    
    ''' Panel for holding the configuration options and canvas panels
    '''
    
    def __init__(self, parent):
        
        ''' Set up the panel, add the controls, bind the events
        '''
        
        casie.gui.panels.base.Panel.__init__(self, parent)
        self.active_sketch = None
    
    def add_controls(self):
        
        ''' Add wx controls to the panel
        '''  
        
        self.panel_label = wx.StaticText(parent = self, label="Sketcher")
        
        self.point_button = wx.Button(parent = self, label="Point")
        self.line_button = wx.Button(parent = self, label="Line")
        self.rectangle_button = wx.Button(parent = self, label="Rectangle")
        self.circle_button = wx.Button(parent = self, label="Circle")
        self.arc_button = wx.Button(parent = self, label="Arc")
        
        self.grid_label = wx.StaticText(parent = self, label="Grid")
        self.grid_on_button = wx.Button(parent = self, label="On")
        self.grid_off_button = wx.Button(parent = self, label="Off")
        
        self.snap_label = wx.StaticText(parent = self, label="Snap")
        self.snap_on_button = wx.Button(parent = self, label="On")
        self.snap_off_button = wx.Button(parent = self, label="Off")
        
        self.undo_button = wx.Button(parent = self, label="Undo")
        self.clear_button = wx.Button(parent = self, label="Clear")
        self.finish_button = wx.Button(parent = self, label="Finish")
        
    def size_controls(self):
        
        ''' Size the controls
        ''' 
        
        expand_fmt = [1, wx.EXPAND, 3]
        centre_fmt = [1, wx.CENTER, 3]
        
        buttonSizer = wx.FlexGridSizer(5, cols=5, vgap=3, hgap=3)
        buttonSizer.Add((0,0), *expand_fmt)  
        buttonSizer.Add((0,0), *expand_fmt)  
        buttonSizer.Add(self.panel_label, *centre_fmt) 
        buttonSizer.Add((0,0), *expand_fmt) 
        buttonSizer.Add((0,0), *expand_fmt)  
        buttonSizer.Add(self.point_button, *expand_fmt)  
        buttonSizer.Add(self.line_button, *expand_fmt) 
        buttonSizer.Add(self.rectangle_button, *expand_fmt)  
        buttonSizer.Add(self.circle_button, *expand_fmt) 
        buttonSizer.Add(self.arc_button, *expand_fmt)
        buttonSizer.Add(self.grid_label, *centre_fmt)  
        buttonSizer.Add(self.grid_on_button, *expand_fmt) 
        buttonSizer.Add(self.grid_off_button, *expand_fmt)
        buttonSizer.Add((0,0), *expand_fmt)  
        buttonSizer.Add((0,0), *expand_fmt)  
        buttonSizer.Add(self.snap_label, *centre_fmt) 
        buttonSizer.Add(self.snap_on_button, *expand_fmt) 
        buttonSizer.Add(self.snap_off_button, *expand_fmt)
        buttonSizer.Add((0,0), *expand_fmt)  
        buttonSizer.Add((0,0), *expand_fmt)  
        buttonSizer.Add(self.undo_button, *expand_fmt)
        buttonSizer.Add(self.clear_button, *expand_fmt)
        buttonSizer.Add(self.finish_button, *expand_fmt)
        
        hbox = wx.BoxSizer()
        hbox.AddSpacer(10)
        hbox.Add(buttonSizer, *expand_fmt)
        hbox.AddSpacer(10)
        
        vbox = wx.BoxSizer(wx.VERTICAL)
        vbox.AddSpacer(10)
        vbox.Add(hbox, 1, wx.CENTER)
        vbox.AddSpacer(10)

        self.SetSizer(vbox)
        vbox.Fit(self)   
        
    def bind_events(self):
        
        ''' Bind the events
        ''' 
        
        self.point_button.Bind(wx.EVT_BUTTON, self.on_draw_point)
        self.line_button.Bind(wx.EVT_BUTTON, self.on_draw_line)
        self.undo_button.Bind(wx.EVT_BUTTON, self.on_undo)
        self.clear_button.Bind(wx.EVT_BUTTON, self.on_clear)
        self.finish_button.Bind(wx.EVT_BUTTON, self.on_finish)
        
    def set_active_sketch(self, sketch):
        
        ''' Set the active sketch
        ''' 
        
        self.active_sketch = sketch
        
    def on_draw_point(self, event):
        
        ''' Set draw mode to point
        '''   
        
        if self.active_sketch is not None:
            self.active_sketch.set_all_false()
            self.active_sketch.draw_point = True
            
    def on_draw_line(self, event):

        ''' Set draw mode to line
        '''   
        
        if self.active_sketch is not None:
            self.active_sketch.set_all_false()
            self.active_sketch.draw_line = True
    
    def on_finish(self, event):
        
        ''' Finish the existing sketch
        '''   
        
        if self.active_sketch is not None:
            if self.active_sketch.type == "Base":
                self.active_sketch.clear() 
                self.refresh_scene_2d_canvas()
            else:
                self.get_canvas_manager().finish_sketch()         
    
    def on_clear(self, event):
        
        ''' Clear the existing sketch
        '''  
        
        if self.active_sketch is not None:
            self.active_sketch.clear()
            self.refresh_scene_2d_canvas()

    def on_undo(self, event):
        
        ''' Undo the last operation in the sketch
        '''  
        
        if self.active_sketch is not None:
            self.active_sketch.undo()
            self.refresh_scene_2d_canvas()
