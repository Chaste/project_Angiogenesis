"""
Test conversion of vessel networks between difference geometrical representations
and storage formats
"""

import unittest
import chaste.image.image_to_surface
import chaste.utility.readwrite
import chaste.simulation.setup

class TestImageToSurface(unittest.TestCase):
    
    def test_bio_2d(self):
        file_handler = chaste.simulation.setup.setup("/home/grogan/test/image/TestImageToSurface/")
#        image = chaste.utility.readwrite.read("/home/grogan/Chaste-Projects/Angiogenesis/test/python/data/bio.tif")
        image = chaste.utility.readwrite.read("/home/grogan/proc_1.tif")
        
        tool = chaste.image.image_to_surface.VtkImageToPolyData2d()
        tool.input = image
        tool.update()
        
        chaste.utility.readwrite.write(image, file_handler.GetOutputDirectoryFullPath() + "/original2d_large.vti")
        chaste.utility.readwrite.write(tool.output, file_handler.GetOutputDirectoryFullPath() + "/vtk_surface2d_large.vtp")
        
    def test_bio_3d(self):
        
        pass
        
if __name__ == '__main__':
    unittest.main()