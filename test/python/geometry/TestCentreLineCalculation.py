import unittest
import chaste.geometry
import chaste.geometry.centreline_calculators2d
import chaste.geometry.converters as converters
import chaste.simulation.setup
import chaste.utility.readwrite

class TestCentreLines2d(unittest.TestCase):
    
    def test_run(self):
        file_handler = chaste.simulation.setup.setup("/home/grogan/test/geometry/TestCentreLines2d/")
        boundary = chaste.utility.readwrite.read("/home/grogan/Chaste-Projects/Angiogenesis/test/python/data/bio_boundaries_extended.vtp")
        
        tool = chaste.geometry.boundary_extractors.BoundaryExtractor2d()
        tool.labels = [1, 2]
        tool.input = boundary
        tool.update()
        edges = tool.output
        
        tool2 = chaste.geometry.surface_extensions.SurfaceExtension2d()
        tool2.input = boundary
        tool2.edges = edges
        tool2.update()
        chaste.utility.readwrite.write(tool2.output, file_handler.GetOutputDirectoryFullPath() + "/bio_boundaries_extended.vtp")
        # 
if __name__ == '__main__':
    unittest.main()