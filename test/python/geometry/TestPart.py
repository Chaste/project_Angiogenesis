''' Tests for the geometry module.
'''

import math
import numpy as np
import unittest
import chaste.geometry

class TestPart(unittest.TestCase):
    
    ''' Test Part Functionality'''
          
    def test_all_methods(self):
        
        # Make a composite Part, a circle in a square
        part = chaste.geometry.Part()
        part.AddRectangle(1.0, 1.0)
        centre = (0.5, 0.5 ,0.0)
        part.AddCircle(0.33, centre)
        
        # Get the VTK Representation
        vtk_part = part.GetVtk(True)