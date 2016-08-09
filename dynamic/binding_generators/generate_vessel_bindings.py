#!/usr/bin/env python

"""
This scipt automatically generates Python bindings using a rule based approach
"""
import os, sys
from pyplusplus import module_builder
from pyplusplus.module_builder import call_policies
from pygccxml import parser

def generate_wrapper_files(args):
    module_name = args[1]
    work_dir = args[2]
    header_collection = args[3]
    includes = args[4:]
    
    generator_path = "/home/grogan/Downloads/castxml/bin/castxml"
    generator_name = "castxml"
    xml_generator_config = parser.xml_generator_configuration_t(xml_generator_path=generator_path, 
                                                                xml_generator=generator_name,
                                                                compiler = "gnu",
                                                                compiler_path="/usr/bin/gcc",
                                                                include_paths=includes)
     
    builder = module_builder.module_builder_t([header_collection],
                                                xml_generator_path = generator_path,
                                                xml_generator_config = xml_generator_config,
                                                include_paths = includes)

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

    builder.global_ns.namespace('std').exclude()
    builder.build_code_creator(module_name="_chaste_project_PyChaste_" + module_name)
    builder.code_creator.user_defined_directories.append( os.path.abspath('.') )
    builder.write_module(work_dir + "/dynamic/" + module_name + "_auto.cpp")

if __name__=="__main__":
    generate_wrapper_files(sys.argv)