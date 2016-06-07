''' Tests for the geometry module.
'''

import unittest
import math
import numpy as np
import chaste.geometry
        
class TestFacet(unittest.TestCase):
    
    ''' Test Facet Functionality'''
          
    def test_all_methods(self):
        
        # Make some vertices
        vertex1 = chaste.geometry.Vertex((0.0, 0.0, 0.0))
        vertex2 = chaste.geometry.Vertex((1.0, 0.0, 0.0))
        vertex3 = chaste.geometry.Vertex((1.0, 1.0, 0.0))
        vertex4 = chaste.geometry.Vertex((0.0, 1.0, 0.0))
    
        # Make a polygon with several vertices
        polygon = chaste.geometry.Polygon([vertex1, vertex2, vertex3, vertex4])

        # Make a facet
        facet = chaste.geometry.Facet(polygon)
        self.assertEqual(len(facet.GetVertices()), 4)
        self.assertEqual(len(facet.GetPolygons()), 1)

        # Check the geometric features
        centroid = (0.5, 0.5, 0.0)
        normal = (0.0, 0.0, 1.0)
            
        # Check translating and rotating
        translation_vector = (1.0, 2.0, 3.0)
        facet.Translate(translation_vector)
        rotation_axis = (0.0, 0.0, 1.0)
        facet.RotateAboutAxis(rotation_axis, math.pi)