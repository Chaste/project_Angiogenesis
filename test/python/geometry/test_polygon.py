''' Tests for the geometry module.
'''

import math
import numpy as np
from nose.tools import assert_equals, assert_almost_equals
from casie.utility import assert_almost_equals_lists
import casie.geometry 
        
class TestPolygon():
    
    ''' Test Polygon Functionality'''
    
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
        
        # Make a polygon with one vertex
        polygon1 = casie.geometry.Polygon(vertex1)
        assert_equals(len(polygon1.get_vertices()), 1)
        polygon1.add_vertices([vertex2, vertex3, vertex4])
        
        # Make a polygon with several vertices
        polygon2 = casie.geometry.Polygon([vertex2, vertex3, vertex4])
        assert_equals(len(polygon2.get_vertices()), 3)
        polygon2.add_vertex(vertex1)
        
        # Check the geometric features
        centroid = (0.5, 0.5, 0.0)
        normal = (0.0, 0.0, 1.0)
        assert_almost_equals_lists(polygon2.get_centroid(), centroid)
        assert_almost_equals_lists(polygon2.get_normal(), normal)
            
        bounding_box = polygon2.get_bbox()
        target_bbox = [0.0, 1.0, 0.0, 1.0, 0.0, 0.0]
        assert_almost_equals_lists(bounding_box, target_bbox)
        
        # Check translating and rotating
        translation_vector = (1.0, 2.0, 3.0)
        polygon2.translate(translation_vector)
        rotation_axis = (0.0, 0.0, 1.0)
        polygon2.rotate(rotation_axis, math.pi)
        