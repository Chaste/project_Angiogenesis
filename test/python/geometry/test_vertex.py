''' Tests for the geometry module.
'''

import math
import numpy as np
from nose.tools import assert_equals, assert_almost_equals
from casie.utility import assert_almost_equals_lists
import casie.geometry

class TestVertex():
    
    ''' Test Vertex Functionality'''
    
    @classmethod
    def setup_class(cls):
        pass
        
    @classmethod
    def teardown_class(self):
        pass
          
    def test_all_methods(self):
        
        # Make a vertex at the specified location
        input_location = np.array((0.0, 1.0, 2.0))
        vertex = casie.geometry.Vertex(input_location)
        
        # Assert that the location is returned correctly
        assert_almost_equals_lists(input_location, vertex.get_location())
        
        # Set the Id and check it
        vertex.id = 10
        assert_equals(vertex.id, 10)
        
        # Move the vertex and check the new location
        translation_vector = np.array((1.0, 2.0, 3.0))
        vertex.translate(translation_vector)
        assert_almost_equals_lists(translation_vector + input_location, vertex.get_location())  
            
        # Rotate the vertex and check the location
        rotation_axis = (0.0, 0.0, 1.0)
        vertex.rotate(rotation_axis, math.pi)
        assert_almost_equals(-(input_location[0] + translation_vector[0]), vertex.get_location()[0], 6)    
        assert_almost_equals(-(input_location[1] + translation_vector[1]), vertex.get_location()[1], 6)    
        assert_almost_equals(input_location[2] + translation_vector[2], vertex.get_location()[2], 6)  