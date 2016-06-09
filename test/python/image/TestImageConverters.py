"""
Test conversion of vessel networks between difference geometrical representations
and storage formats
"""

import unittest
import numpy as np
import chaste.image.converters as converters
import chaste.mesh
import chaste.simulation.setup
import chaste.utility.readwrite
import scipy.misc

class TestVtkImageDataToNumpy(unittest.TestCase):
    
    def test_run(self):
        file_handler = chaste.simulation.setup.setup("/home/grogan/test/TestVtkImageDataToNumpy/")
        
        # Set up a vtk solution
        grid = chaste.mesh.RegularGrid()
        grid.SetSpacing(1.0)
        grid.SetExtents((20, 20, 20))
        centre = np.array((10.0, 10.0, 10.0))
        
        values = []
        for idx in range(grid.GetNumberOfPoints()): 
            eachLocation = grid.GetLocationOf1dIndex(idx)
            values.append(np.linalg.norm(np.array(eachLocation-centre)))
        
        grid.SetUpVtkGrid()
        grid.SetPointValues(values)
        grid.Write(file_handler)
        
        grid = chaste.utility.readwrite.read(file_handler.GetOutputDirectoryFullPath() + "/grid.vti")
        converter = converters.VtkImageDataToNumpy()
        converter.input = grid
        converter.update()
        numpy_array = converter.output
        
        scipy.misc.imsave(file_handler.GetOutputDirectoryFullPath() + "/image.jpg", numpy_array)
        # 
if __name__ == '__main__':
    unittest.main()