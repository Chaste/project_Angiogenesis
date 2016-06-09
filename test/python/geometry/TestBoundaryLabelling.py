import unittest
import chaste.geometry
import chaste.geometry.boundary_markers
import chaste.geometry.converters as converters
import chaste.simulation.setup
import chaste.utility.readwrite

class TestBoundaryLabelling2d(unittest.TestCase):
    
    def test_run(self):
        file_handler = chaste.simulation.setup.setup("/home/grogan/test/geometry/TestBoundaryLabelling2d/")
        boundary = chaste.utility.readwrite.read("/home/grogan/Chaste-Projects/Angiogenesis/test/python/data/bio_boundaries.vtp")
        
        tool = chaste.geometry.boundary_markers.BoundaryMarker2d()
        tool.seed_points = [[0, 95.0, 0.0], [0.0, 137.0, 0.0], [22.0, 163.0, 0.0], [50.0, 163.0, 0.0], [65.0, 0.0, 0.0], [130.0, 0.0, 0.0]]
        tool.labels = [1, 1, 1, 1, 2, 2]
        tool.input = boundary
        tool.update()
        
        chaste.utility.readwrite.write(tool.output, file_handler.GetOutputDirectoryFullPath() + "/labelled_boundary.vtp")
        # 
if __name__ == '__main__':
    unittest.main()