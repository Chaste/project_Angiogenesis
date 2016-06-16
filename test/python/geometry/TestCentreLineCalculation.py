import unittest
import chaste.geometry
import chaste.geometry.centreline_calculators2d
import chaste.geometry.converters as converters
import chaste.simulation.setup
import chaste.utility.readwrite

class TestCentreLines2d(unittest.TestCase):
    
    def test_run(self):
        file_handler = chaste.simulation.setup.setup("/home/grogan/test/geometry/TestCentreLines2d/")
#        boundary = chaste.utility.readwrite.read("/home/grogan/Chaste-Projects/Angiogenesis/test/python/data/bio_boundaries_extended.vtp")
        
        boundary = chaste.utility.readwrite.read("/home/grogan/test/image/TestImageToSurface/vtk_surface2d.vtp")
        
        tool = chaste.geometry.centreline_calculators2d.Centrelines2d()
        tool.surface = boundary
        tool.update()

        chaste.utility.readwrite.write(tool.network, file_handler.GetOutputDirectoryFullPath() + "/centreline.vtp")
        # 
if __name__ == '__main__':
    unittest.main()