"""
Show a tiff using vtk
"""

import code.image.settings
from code.image.tools import simple_render
from code.image.tools import it_wrapper

work_dir, input_data = code.image.settings.do_setup("show")
input_path = work_dir + "/threshold/itk_tile_0_2_6.tiff"

image_type = it_wrapper.default_image_type()
reader = it_wrapper.read_image(input_path, image_type)
correct = it_wrapper.correct_spacing(reader, image_type, input_data["VoxelSizeZ"])
converter = it_wrapper.convert_vtk(reader, image_type)
    
simple_render.render_image(converter, threshold = 1.0)