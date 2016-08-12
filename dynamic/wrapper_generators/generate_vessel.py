#!/usr/bin/env python

"""
This scipt automatically generates Python bindings using a rule based approach
"""
import sys
from pyplusplus import module_builder
from pyplusplus.module_builder import call_policies
from pygccxml import parser

def update_builder(builder):

#    builder.print_declarations()
#     include_classes = ["NodeFlowProperties< 3 >", 
#                        "SegmentFlowProperties< 3 >", 
#                        "VesselFlowProperties< 3 >"]
    include_classes = ["NodeFlowProperties< 3 >"]#, "SegmentFlowProperties< 3 >"]
    for eachClass in include_classes:
        builder.class_(eachClass).include()
        
    builder.class_('NodeFlowProperties< 3 >').rename('NodeFlowProperties3')
    
#    builder.class_('SegmentFlowProperties< 3 >').rename('SegmentFlowProperties')
#    builder.class_('VesselFlowProperties< 3 >').rename('VesselFlowProperties')
        
#    builder.member_functions('OpenOutputFile').exclude()
#    builder.member_functions('FindMatches').exclude()

    return builder