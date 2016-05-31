import casie.geometry
import casie.plot.three.scene
reload(casie.geometry)
reload(casie.plot.three)

# Create a circle with radius 10
part = casie.geometry.Part()
circle = part.AddCircle(100)
part.Extrude(circle, 100)

# Visualize it using VTK
scene = casie.plot.three.scene.Scene()
scene.add_part(part)
scene.show(interactive=True)