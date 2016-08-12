#!/usr/bin/env python

"""
This scipt automatically generates Python bindings using a rule based approach
"""
import sys
from pyplusplus import module_builder
from pyplusplus.module_builder import call_policies
from pygccxml import parser
import generate_flow
import generate_utility
import generate_vessel

def generate_wrappers(args):
    module_name = args[1]
    work_dir = args[2]
    header_collection = args[3]
    castxml_binary = args[4]
    includes = args[5:]
    
    xml_generator_config = parser.xml_generator_configuration_t(xml_generator_path=castxml_binary, 
                                                                xml_generator="castxml",
                                                                compiler = "gnu",
                                                                compiler_path="/usr/bin/gcc",
                                                                include_paths=includes)
     
    builder = module_builder.module_builder_t([header_collection],
                                                xml_generator_path = castxml_binary,
                                                xml_generator_config = xml_generator_config,
                                                include_paths = includes)
    
    if("flow" in module_name):
        builder = generate_flow.update_builder(builder)
        
    if("utility" in module_name):
        builder = generate_utility.update_builder(builder)
        
    if("vessel" in module_name):
        builder = generate_vessel.update_builder(builder)

    builder.global_ns.namespace('std').exclude()
    builder.build_code_creator(module_name="_chaste_project_PyChaste_" + module_name)
    builder.write_module(work_dir + "/dynamic/" + module_name + ".cpp")

if __name__=="__main__":
    generate_wrappers(sys.argv)
