""" 
Properties for the chaste GUI.
"""

try:
    import wx
    _have_wx = True
except ImportError:
    _have_wx = False
    
_colors = {"background": (0, 0, 0, 1.0),
           "panel_background": (255, 255, 255, 1.0),
          "none": (192, 192, 192, 0.75)}

_tools = [("ModelBuilder", "chaste.gui.model_builder.ModelBuilder"),
          ("BoundaryMarker2d", "chaste.geometry.gui.labelling.BoundaryMarker2d"), 
          ("Centrelines2d", "chaste.geometry.gui.centrelines.Centrelines2d"),
          ("SurfaceFromImage2d", "chaste.image.converter.SurfaceFromImage2d"),
          ("ImageViewer2d", "chaste.image.viewer.ImageViewer2d"),
          ("Mesh2d", "chaste.mesh.meshers.Mesh2d")]