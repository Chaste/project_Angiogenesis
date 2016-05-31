''' The Casie GUI app and main frame.

'''

import wx
import vtk

import casie.gui.sub_panels.scene3d_canvas
import casie.geometry

class Frame(wx.Frame):
    
    def __init__(self, parent):
        
        '''
        Set up the main frame.
        
        Create managers for event handling and the main menus and panels.
        ''' 
        
        wx.Frame.__init__(self, parent, size = (800,600))

        self.create_panels()
        self.Show(True)
    
    def create_panels(self):
        
        ''' 
        Create and size the main panels
        ''' 
        
        self.panel = wx.Panel(self, wx.ID_ANY)
        self.panel.choice = casie.gui.sub_panels.scene3d_canvas.Panel(self.panel)
        
        #domain = casie.geometry.Part()
        #domain.AddCuboid(100.0, 100.0, 1.0)
        
        #polyd = domain.GetVtk()
       # mapper = vtk.vtkPolyDataMapper()
        #mapper.SetInput(polyd)
        #mapper.ScalarVisibilityOn()
        #actor = vtk.vtkActor()
        #actor.SetMapper(mapper)
        
        #self.panel.choice.add_mesh(actor)

        hbox = wx.BoxSizer()
        hbox.Add(self.panel.choice, 1, wx.EXPAND, 1)
        self.panel.SetSizer(hbox) 
                
if __name__ == "__main__":

    # Launch the window
    app = wx.App(redirect=False)
    app.frame = Frame(parent=None)
    app.MainLoop()