''' 3D Vtk Graphics Window
'''

import wx
import vtk
from vtk.wx.wxVTKRenderWindowInteractor import wxVTKRenderWindowInteractor

class Panel(wx.Panel):
    
    ''' Vtk panel for 3D plotting
    '''
    
    def __init__(self, parent):
        wx.Panel.__init__(self, parent)
         
        self.widget = wxVTKRenderWindowInteractor(self, -1)
        self.widget.Enable(1)
        self.widget.AddObserver("ExitEvent", lambda o,e,f=self: f.Close())
        
        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.sizer.Add(self.widget, 1, wx.LEFT | wx.TOP | wx.EXPAND)
        self.SetSizer(self.sizer)
        self.Layout()
        
        self.initialize_display()
        
    def initialize_display(self):  
        
        self.ren = vtk.vtkRenderer()
        self.ren.SetBackground(1.0, 1.0, 1.0)
        self.widget.GetRenderWindow().AddRenderer(self.ren)
        self.widget.GetRenderWindow().SetAlphaBitPlanes(1)
                
        axes = vtk.vtkAxesActor()
        self.marker = vtk.vtkOrientationMarkerWidget()
        self.marker.SetInteractor(self.widget._Iren )
        self.marker.SetOrientationMarker( axes )
        self.marker.SetViewport(0.75,0,1,0.25)
        self.marker.SetEnabled(1)
        
        self.textActor = vtk.vtkTextActor()
        self.textActor.SetTextScaleModeToProp()
        self.textActor.SetDisplayPosition(50, 25)
        
        self.reset_camera()
        
    def reset_camera(self):
        
        self.ren.ResetCamera()
        self.ren.ResetCameraClippingRange()
        cam = self.ren.GetActiveCamera()
        cam.Elevation(10)
        cam.Azimuth(70) 
        
    def add_cell_population(self, actor):
        
        self.ren.AddActor(actor) 
        self.reset_camera()
        
    def add_vessel_network(self, actor):
        
        self.ren.AddActor(actor[0]) 
        self.ren.AddActor(actor[1]) 
        self.reset_camera()    
        
    def add_part(self, actor):
         
        self.ren.AddActor(actor[0]) 
        self.ren.AddActor(actor[1]) 
        self.reset_camera()
        
    def add_mesh(self, actor):
        
#         self.ren.AddActor(actor[0]) 
#         self.ren.AddActor(actor[1]) 
        self.ren.AddActor(actor) 
        self.reset_camera()     
        
    def add_image(self, image):
         
        threshold = 80.0
        alphaChannelFunc = vtk.vtkPiecewiseFunction()
        alphaChannelFunc.AddPoint(0, 0.0)
        alphaChannelFunc.AddPoint(threshold, 0.0)
        alphaChannelFunc.AddPoint(threshold + 1, 1.0)
        alphaChannelFunc.AddPoint(255, 1.0)
         
        colorFunc = vtk.vtkColorTransferFunction()
        colorFunc.AddRGBPoint(0, 0.0, 0.0, 1.0)
        colorFunc.AddRGBPoint(threshold, 1.0, 0.0, 0.0)
        colorFunc.AddRGBPoint(255, 1.0, 0.0, 0.0)
        volumeProperty = vtk.vtkVolumeProperty()
        volumeProperty.SetColor(colorFunc)
        volumeProperty.SetScalarOpacity(alphaChannelFunc)
           
        compositeFunction = vtk.vtkVolumeRayCastCompositeFunction()
        volumeMapper = vtk.vtkVolumeRayCastMapper()
        volumeMapper.SetVolumeRayCastFunction(compositeFunction)
        volumeMapper.SetInputConnection(image.GetProducerPort())
        volume = vtk.vtkVolume()
        volume.SetMapper(volumeMapper)
        volume.SetProperty(volumeProperty)
         
        self.ren.AddVolume(volume)
        
    def save(self, fileName):
        if fileName is not None:
            screenFilter = vtk.vtkWindowToImageFilter()
            screenFilter.SetInput(self.widget.GetRenderWindow())
            screenFilter.SetMagnification(1)
            screenFilter.SetInputBufferTypeToRGBA()
            screenFilter.ReadFrontBufferOff()  
            screenFilter.Update()    
            writer = vtk.vtkPNGWriter()
            writer.SetFileName(fileName + ".png")
            writer.SetInputConnection(screenFilter.GetOutputPort());
            writer.Write() 