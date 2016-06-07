"""
Median filter on chosen image
"""

import logging
from code.image.tools import it_wrapper
import code.image.settings

if __name__ == "__main__":
    
    # Do setup
    tool_name = "median"
    work_dir, image_data = code.image.settings.do_setup(tool_name)
    logger1 = logging.getLogger('processing.'+tool_name)

    output_path = work_dir+"/median"
    
    logger1.info('Start Filtering Image at: ' + work_dir)
    image_type = it_wrapper.default_image_type()
    
    input_path = work_dir + "/convert/itk_tile_0_2_6.tiff"
    logger1.info("Reading the image: " +input_path + " into ITK")
    reader = it_wrapper.read_image(input_path, image_type)
    correct = it_wrapper.correct_spacing(reader, image_type, image_data["VoxelSizeZ"])
            
    logger1.info('Doing Median Filter')
    median_filter = it_wrapper.median_filter(correct, image_type, 1)
    
    logger1.info('Writing TIFF Stack')
    it_wrapper.write_image(output_path+"/itk_tile_0_2_6.tiff", image_type, median_filter)