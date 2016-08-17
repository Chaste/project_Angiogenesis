#!/usr/bin/env python

"""
This scipt automatically generates Python bindings using a rule based approach
"""
import sys
from pyplusplus import module_builder
from pyplusplus.module_builder import call_policies
from pygccxml import parser

def update_builder(builder):

    include_classes = ["NodeFlowProperties< 3 >", 
                       "SegmentFlowProperties< 3 >", 
                       "VesselFlowProperties< 3 >", 
                       "VesselNode< 3 >",
                       "VesselSegment< 3 >",
                       "Vessel< 3 >", 
                       "VesselNetwork< 3 >",
                       "VasculatureGenerator< 3 >",
                       "AbstractVesselNetworkComponent< 3 >",
                       "VtkVesselNetworkWriter< 3 >"]
    for eachClass in include_classes:
        builder.class_(eachClass).include()
        
    builder.class_('NodeFlowProperties< 3 >').rename('NodeFlowProperties3')
    builder.class_('SegmentFlowProperties< 3 >').rename('SegmentFlowProperties3')
    builder.class_('VesselFlowProperties< 3 >').rename('VesselFlowProperties3')
    builder.class_('VesselNode< 3 >').rename('VesselNode3')
    builder.class_('VesselSegment< 3 >').rename('VesselSegment3')
    builder.class_('Vessel< 3 >').rename('Vessel3')
    builder.class_('VesselNetwork< 3 >').rename('VesselNetwork3')      
    builder.class_('VasculatureGenerator< 3 >').rename('VasculatureGenerator3')      
    builder.class_('AbstractVesselNetworkComponent< 3 >').rename('AbstractVesselNetworkComponent3')      
    builder.class_('VtkVesselNetworkWriter< 3 >').rename('VtkVesselNetworkWriter3')     
       

    return builder