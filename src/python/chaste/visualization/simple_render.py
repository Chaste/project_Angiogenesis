"""
Rendering of 3D TIFF stacks with VTK
"""

import vtk

def render_image(image_container, threshold = 200.0):
    
    """
    Render the 3D stack with VTK
    """
    
    volume_mapper = vtk.vtkVolumeRayCastMapper()
    volume_mapper.SetInput(image_container.GetOutput() )
    composite_function = vtk.vtkVolumeRayCastCompositeFunction()
    volume_mapper.SetVolumeRayCastFunction( composite_function )
    
    color_transfer_func = vtk.vtkColorTransferFunction()
    color_transfer_func.AddRGBPoint( 0, 0.0, 1.0, 0.0 )
    color_transfer_func.AddRGBPoint( threshold-1, 0.0, 1.0, 0.0 )
    color_transfer_func.AddRGBPoint( threshold, 1.0, 0.0, 0.0 )
    color_transfer_func.AddRGBPoint( 255.0, 1.0, 0.0, 0.0 )
    
    opacity_transfer_func = vtk.vtkPiecewiseFunction()
    opacity_transfer_func.AddPoint( 0, 0.0 )
    opacity_transfer_func.AddPoint( threshold-1.0, 0.0 )
    opacity_transfer_func.AddPoint( threshold, 1.0 )
    opacity_transfer_func.AddPoint( 255.0, 1.0 )
    
    volume_properties = vtk.vtkVolumeProperty()
    volume_properties.SetColor( color_transfer_func )
    volume_properties.SetScalarOpacity( opacity_transfer_func )
    
    volume = vtk.vtkVolume()
    volume.SetMapper( volume_mapper )
    volume.SetProperty( volume_properties )
    
    renderer = vtk.vtkRenderer()
    renderer.AddVolume( volume )
    
    renderer.ResetCamera()
    
    render_window = vtk.vtkRenderWindow()
    window_interactor = vtk.vtkRenderWindowInteractor()
    
    render_window.AddRenderer( renderer )
    window_interactor.SetRenderWindow( render_window )

    render_window.Render()
    window_interactor.Start()