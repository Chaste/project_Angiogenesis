"""
Convert to a vtk image
"""

import logging
from vmtk import pypes
import code.image.settings

if __name__ == "__main__":
    
    # Do setup
    tool_name = "vtk_image"
    work_dir, image_data = code.image.settings.do_setup(tool_name)
    logger1 = logging.getLogger('processing.'+tool_name)

    output_path = work_dir+"/vtk_image"
    
    logger1.info('Start VTK Converting Image at: ' + work_dir)

    #input_path = work_dir + "/threshold/itk_tile_0_2_6.tiff"
    input_path = "/home/grogan/clipped.tif"
    output_path = "/home/grogan/clipped.vti"
    myArguments = 'vmtkimagereader -ifile ' + input_path + ' --pipe vmtkimagewriter -ofile ' + output_path

    myPype = pypes.PypeRun(myArguments)