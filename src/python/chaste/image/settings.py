"""
Settings for image processing scripts
"""

import os
import logging
from code.image.processing.conversion import bftools_wrapper

def do_setup(tool_name = None):
    
    work_dir = os.environ['IMAGE_WORK_DIR'] + "/NormalVasculature/"
    
    if tool_name is not None:
        filename = work_dir + "/logging/" + tool_name + ".log"
    else:
        filename = work_dir + "/logging/root.log"
        
    logging.basicConfig(level=logging.DEBUG,
                format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                datefmt='%m-%d %H:%M',
                filename=filename,
                filemode='w')
    
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    
    # set a format which is simpler for console use
    formatter = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')
    
    # tell the handler to use this format
    console.setFormatter(formatter)
    
    # add the handler to the root logger
    logging.getLogger('').addHandler(console)
    
    # try to get the image metadata
    try:
        f = open(work_dir + "/metadata/raw.dat", "r")
    except:
        logging.warning("Image metadata file not found. Attempting to generate one.")
        bftools_wrapper.get_info(work_dir+"/raw/image.lsm", work_dir+"/metadata/raw.dat")
        f = open(work_dir + "/metadata/raw.dat", "r")
        
    image_data = read_metadata(f)
    f.close()
        
    return [work_dir, image_data]

def read_metadata(f):
    
    image_data = {}
    image_data["VoxelSizeX"] = -1.0
    image_data["VoxelSizeY"] = -1.0
    image_data["VoxelSizeZ"] = -1.0
    vox_x_keyword = 'VoxelSizeX' 
    vox_y_keyword = 'VoxelSizeY' 
    vox_z_keyword = 'VoxelSizeZ' 
        
    image_data["DimensionX"] = -1.0
    image_data["DimensionY"] = -1.0
    image_data["DimensionZ"] = -1.0
    dim_x_keyword = 'DimensionX' 
    dim_y_keyword = 'DimensionY' 
    dim_z_keyword = 'DimensionZ'     
        
    for line in f:
        if line.startswith(vox_x_keyword):
            image_data["VoxelSizeX"]  = float(line.split(':')[1])
        if line.startswith(vox_y_keyword):
            image_data["VoxelSizeY"] = float(line.split(':')[1])
        if line.startswith(vox_z_keyword):
            image_data["VoxelSizeZ"] = float(line.split(':')[1])
        if line.startswith(dim_x_keyword):
            image_data["DimensionX"] = int(line.split(':')[1])
        if line.startswith(dim_y_keyword):
            image_data["DimensionY"] = int(line.split(':')[1])
        if line.startswith(dim_z_keyword):
            image_data["DimensionZ"] = int(line.split(':')[1])
    
    return image_data