#!/usr/bin/env python

"""
This scipt automatically generates Python bindings using a rule based approach
"""
import sys
from pyplusplus import module_builder
from pyplusplus.module_builder import call_policies
from pygccxml import parser
import generate_bindings

def update_builder(builder):

    include_classes = ["NodeFlowProperties<3>", 
                       "SegmentFlowProperties<3>", 
                       "VesselFlowProperties<3>", 
                       "VesselNode<3>",
                       "VesselSegment<3>",
                       "Vessel<3>", 
                       "VesselNetwork<3>",
                       "VasculatureGenerator<3>",
                       "AbstractVesselNetworkComponent<3>",
                       "VtkVesselNetworkWriter<3>",
                       "VtkVesselNetworkReader<3>",
                       "VesselDistribution",
                       "SegmentLocation"]
    
    for eachClass in include_classes:
        builder.class_(eachClass).include()  
        new_name = generate_bindings.template_replace(eachClass)
        if(new_name != eachClass):
            builder.class_(eachClass).rename(new_name) 
            
    # Default template args problem
    builder.class_('VasculatureGenerator<3> ').calldefs().use_default_arguments=False    
    
    # Custom constructors
    builder.add_declaration_code('boost::shared_ptr<VesselNode<3> > (*VN3_Doubles)(double, double, double) = &VesselNode<3>::Create;')
 #   builder.add_declaration_code('boost::shared_ptr<VesselNode<3> > (*VN3_Doubles)(double, double, double) = &VesselNode<3>::Create;')
    builder.class_('VesselNode<3>').add_registration_code('def("__init__", bp::make_constructor(VN3_Doubles))')


    return builder