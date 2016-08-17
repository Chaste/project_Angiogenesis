#!/usr/bin/env python

"""
This scipt automatically generates Python bindings using a rule based approach
"""
import sys
from pyplusplus import module_builder
from pyplusplus.module_builder import call_policies
from pygccxml import parser

def update_builder(builder):

    include_classes = ["RegularGrid<3, 3>", "HybridMesh<3, 3>", "SharedPottsMeshGenerator<3>", "PottsMesh<3>" ]
    for eachClass in include_classes:
        builder.class_(eachClass).include()

    builder.class_('RegularGrid< 3, 3 >').rename('RegularGrid3')
    builder.class_('HybridMesh< 3, 3 >').rename('HybridMesh3')
    builder.class_('HybridMesh< 3, 3 >').member_functions("GenerateFromStl").exclude()
    builder.class_('HybridMesh< 3, 3 >').member_functions("GenerateTriMeshFromPolyData").exclude()
    builder.class_('PottsMesh<3>').member_functions("GetElement").exclude()
   
    return builder