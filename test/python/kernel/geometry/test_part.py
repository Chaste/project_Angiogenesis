''' Tests for the geometry module.
'''

import math
import numpy as np
from nose.tools import assert_equals, assert_almost_equals
from casie.utility import assert_almost_equals_lists
import casie.geometry

class TestPart():
    
    ''' Test Part Functionality'''
    
    @classmethod
    def setup_class(cls):
        pass
        
    @classmethod
    def teardown_class(self):
        pass
          
    def test_all_methods(self):
        
        # Make a composite Part, a circle in a square
        part = casie.geometry.Part()
        part.AddRectangle(1.0, 1.0)
        centre = (0.5, 0.5 ,0.0)
        part.AddCircle(0.33, centre)
        
        # Get the VTK Representation
        vtk_part = part.GetVtk(True)
        print vtk_part
        