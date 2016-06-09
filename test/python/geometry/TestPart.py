''' Tests for the geometry module.
'''

import math
import numpy as np
import unittest
import chaste.geometry
import chaste.core

class TestPart(unittest.TestCase):
    
    ''' Test Part Functionality'''
          
    def test_all_methods(self):
        
        file_handler = chaste.core.OutputFileHandler("geometry/TestPart", True)
        
        # Make a composite Part, a circle in a square
        part = chaste.geometry.Part()
        part.AddRectangle(1.0, 1.0)
        centre = (0.5, 0.5 ,0.0)
        part.AddCircle(0.33, centre)
        
        # Get the VTK Representation
        part.Write(file_handler.GetOutputDirectoryFullPath() + "original_part.vtp")
        
        vtk_part = part.GetVtk(True)
        
if __name__ == '__main__':
    unittest.main()
        
        
        