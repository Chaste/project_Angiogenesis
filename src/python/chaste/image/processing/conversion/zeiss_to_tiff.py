"""
Conversion of LSM/CZI files to OME and Generic TIFF
"""

import logging
from code.image.tools import bftools_wrapper
from code.image.tools import it_wrapper
import code.image.settings
    
if __name__ == "__main__":
    
    # Do setup
    tool_name = "convert"
    work_dir, image_data = code.image.settings.do_setup(tool_name)
    logger1 = logging.getLogger('processing.'+tool_name)

    input_path = work_dir + "/raw/image.lsm"
    output_path = work_dir+"/convert"
    
    logger1.info('Start Converting Image at: ' + input_path)
    
    # Check if the image needs to be tiled
    max_dimension = 4096
    tile_x = -1
    tile_y = -1
    num_tiles_x = -1
    num_tiles_y = -1
    
    if image_data["DimensionX"] > max_dimension:
        tile_x = max_dimension
    if image_data["DimensionY"] > max_dimension:
        tile_y = max_dimension

    if tile_x == max_dimension and tile_y == -1:
        tile_y = max_dimension
        
    if tile_y == max_dimension and tile_x == -1:
        tile_x = max_dimension
        
    if tile_x > -1:
        num_tiles_x = int(float(image_data["DimensionX"])/(max_dimension)) + 1
        num_tiles_y = int(float(image_data["DimensionY"])/(max_dimension)) + 1
            
        logger1.warning('Max image dimensions exceeded, tiling before further processing.')
        bftools_wrapper.convert(input_path, output_path, tile_x, tile_y, clip_range = None, channel=-1)
        logger1.info('Tiling completed.')
    
    # Set the iamge type
    image_type = it_wrapper.default_image_type()

    if tile_x > -1:
        for jdx in range(num_tiles_y):
            for idx in range(num_tiles_x):
                input_path = work_dir + "/convert/ome_tile_" + str(idx) + "_" + str(jdx) +"_"+ str(idx+ num_tiles_x*jdx) + ".tiff"
                logger1.info("Reading the image: " +input_path + " into ITK")
                reader = it_wrapper.read_image(input_path, image_type)
                logger1.info('Writing TIFF Stack')
                correct = it_wrapper.correct_spacing(reader, image_type, image_data["VoxelSizeZ"])
                it_wrapper.write_image(output_path+"/itk_tile_" + str(idx) + "_" + str(jdx) +"_"+ str(idx+ num_tiles_x*jdx) + ".tiff", image_type, correct)
    else:
        # Read the image
        logger1.info('Reading the image into ITK')
        reader = it_wrapper.read_image(input_path + "/convert/ome.tiff", image_type)
        correct = it_wrapper.correct_spacing(reader, image_type, image_data["VoxelSizeZ"])
    
        # Write the image in itk tiff format
        logger1.info('Writing TIFF Stack')
        it_wrapper.write_image(output_path+"/itk.tiff", image_type, correct)
    
    logger1.info('Completed Processing')