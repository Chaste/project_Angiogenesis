"""
Test conversion of vessel networks between difference geometrical representations
and storage formats
"""

import unittest
import chaste.utility.readwrite
import chaste.setup
import chaste.image

class TestImageReader(unittest.TestCase):
    
    def test_bio_2d(self):
        file_handler = chaste.setup.setup("/home/grogan/test/python/image/TestImageReader/")
        reader = chaste.image.ImageReader()
        reader.SetFilename("/home/grogan/Chaste-Projects/Angiogenesis/test/python/data/bio.tif")
        reader.Update()
        chaste.utility.readwrite.write(reader.GetOutput(), file_handler.GetOutputDirectoryFullPath() + "/output_image.vti")
        
class TestImageToSurface(unittest.TestCase):
    
    def test_bio_2d(self):
        file_handler = chaste.setup.setup("/home/grogan/test/python/image/TestImageToSurface/")
        reader = chaste.image.ImageReader()
        reader.SetFilename("/home/grogan/Chaste-Projects/Angiogenesis/test/python/data/bio.tif")
        reader.Update()
        
        surface_extract = chaste.image.ImageToSurface()
        surface_extract.SetInput(reader.GetOutput())
        surface_extract.SetThreshold(1.0, False)
        surface_extract.SetUseMarchingCubes(False)
        surface_extract.Update()
        chaste.utility.readwrite.write(surface_extract.GetOutput(), file_handler.GetOutputDirectoryFullPath() + "/surface_from_image3.vtp")
        
if __name__ == '__main__':
    unittest.main()