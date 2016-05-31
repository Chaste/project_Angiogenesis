''' Tests for the geometry module.
'''

import math
import numpy as np
from nose.tools import assert_equals, assert_almost_equals
from casie.utility import assert_almost_equals_lists
import casie.geometry

        
class TestFacet():
    
    ''' Test Facet Functionality'''
    
    @classmethod
    def setup_class(cls):
        pass
        
    @classmethod
    def teardown_class(self):
        pass
          
    def test_all_methods(self):
        
        # Make some vertices
        vertex1 = casie.geometry.Vertex((0.0, 0.0, 0.0))
        vertex2 = casie.geometry.Vertex((1.0, 0.0, 0.0))
        vertex3 = casie.geometry.Vertex((1.0, 1.0, 0.0))
        vertex4 = casie.geometry.Vertex((0.0, 1.0, 0.0))
    
        # Make a polygon with several vertices
        polygon = casie.geometry.Polygon([vertex1, vertex2, vertex3, vertex4])

        # Make a facet
        facet = casie.geometry.Facet(polygon)
        assert_equals(len(facet.GetVertices()), 4)
        assert_equals(len(facet.GetPolygons()), 1)

        # Check the geometric features
        centroid = (0.5, 0.5, 0.0)
        normal = (0.0, 0.0, 1.0)
        assert_almost_equals_lists(facet.GetCentroid(), centroid)
        assert_almost_equals_lists(facet.GetNormal(), normal)
            
        # Check translating and rotating
        translation_vector = (1.0, 2.0, 3.0)
        facet.Translate(translation_vector)
        rotation_axis = (0.0, 0.0, 1.0)
        facet.RotateAboutAxis(rotation_axis, math.pi)
        