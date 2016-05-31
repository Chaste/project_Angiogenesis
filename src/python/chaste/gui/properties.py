""" 
Properties for the Casie GUI.
"""

try:
    import wx
    _have_wx = True
except ImportError:
    _have_wx = False
    
_colors = {"background": (0, 0, 0, 1.0),
           "panel_background": (255, 255, 255, 1.0),
          "none": (192, 192, 192, 0.75)}

_tools = [("ModelBuilder", "casie.gui.model_builder.ModelBuilder"),
          ("BoundaryMarker2d", "casie.geometry.labelling.BoundaryMarker2d"), 
          ("Centrelines2d", "casie.geometry.centrelines.Centrelines2d"),
          ("SurfaceFromImage2d", "casie.image.converter.SurfaceFromImage2d"),
          ("ImageViewer2d", "casie.image.viewer.ImageViewer2d"),
          ("Mesh2d", "casie.mesh.meshers.Mesh2d")]