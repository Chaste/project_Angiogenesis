"""
Merging of tiled images
"""

import logging
from code.image.tools import it_wrapper
import code.image.settings

if __name__ == "__main__":
    
    # Do setup
    tool_name = "tile"
    work_dir, image_data = code.image.settings.do_setup(tool_name)
    logger1 = logging.getLogger('processing.'+tool_name)

    output_path = work_dir+"/projection"
    
    logger1.info('Start Merging Images at: ' + work_dir)
    image_type = it_wrapper.default_image_type()
    
    num_tiles_y = 3
    num_tiles_x = 3
    container = []
    for jdx in range(num_tiles_y):
        for idx in range(num_tiles_x):
            
            label = str(idx) + "_" + str(jdx) +"_"+ str(idx+ num_tiles_x*jdx)
            input_path = work_dir + "/projection/sum_tile_" + label + ".tiff"
            logger1.info("Reading the image: " +input_path + " into ITK")
            reader = it_wrapper.read_image(input_path, image_type)
            correct = it_wrapper.correct_spacing(reader, image_type, image_data["VoxelSizeZ"])
            container.append(correct)
    
    logger1.info('Merging Tiles')    
    merged = it_wrapper.merge_tiles(container, image_type, num_tiles_x)
      
    output_image_type = it_wrapper.default_image_type()
    logger1.info('Writing TIFF Stack')
    it_wrapper.write_image(output_path+"/sum_merge.tiff", output_image_type, merged)