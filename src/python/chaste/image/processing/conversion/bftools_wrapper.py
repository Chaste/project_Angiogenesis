"""
Functions for tiling and cropping large files in ZEISS formats before further processing
with ITK.
"""

import os
import subprocess

def convert(input_path, output_path, tile_x, tile_y, clip_range = None, channel=-1):
    
    """
    Convert an input image in ZEISS LSM or CZI format to an OME TIFF.
    
    Slice spacing information is lost in this operation.
    
    clip_range crops the image to the specified [x1, y1, x2, y2] dimensions
    
    channel extracts a requested channel
    
    """

    # Set up the command for bfconvert
    command = os.environ['BFTOOLS_DIR'] + "/bfconvert "
    
    log_path = output_path + "/conversion.log"
    output_path += "/ome"
    
    # Clip the image if required
    if clip_range is not None:
        command += " -crop " 
        for idx in range(len(clip_range)-1):
            command += str(clip_range[idx]) + ","
        command += str(clip_range[len(clip_range)-1]) + " "
        output_path += "_cropped_"

    if tile_x > -1:
        command += " -tilex " + str(tile_x) + " -tiley " + str(tile_y) + " "
        output_path += "_tile_%x_%y_%m"
        
    # Extract a specified channel if required
    if channel > -1:
        command += " -channel " + str(channel) + " "
    
    command += input_path + " " + output_path + ".tiff  > " + log_path

    # Do the conversion
    p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
    p.wait()

    
def get_info(input_path, output_path):
    
    """
    Get the image metadata and pipe it to an output text file.
    """
    
    command = os.environ['BFTOOLS_DIR'] +"/showinf -omexml -nopix " + input_path + " > " + output_path
    p = subprocess.Popen(command, shell=True)
    p.wait()
