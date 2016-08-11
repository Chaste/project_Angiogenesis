#!/usr/bin/env python

"""
This scipt automatically generates Python bindings for the utility module.
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
    
    include_classes = ["UnitTester"]#, "SegmentFlowProperties< 3 >"]
    for eachClass in include_classes:
        builder.class_(eachClass).include()

    builder.class_("kg_instance_t<true>").include()
#    builder.class_("boost::units::unit< boost::units::list< boost::units::dim< boost::units::mass_base_dimension, boost::units::static_rational< 1, 1 > >, boost::units::dimensionless_type >, boost::units::homogeneous_system< boost::units::list< boost::units::si::meter_base_unit, boost::units::list< boost::units::scaled_base_unit< boost::units::cgs::gram_base_unit, boost::units::scale< 10, static_rational< 3 > > >, boost::units::list< boost::units::si::second_base_unit, boost::units::list< boost::units::si::ampere_base_unit, boost::units::list< boost::units::si::kelvin_base_unit, boost::units::list< boost::units::si::mole_base_unit, boost::units::list< boost::units::si::candela_base_unit, boost::units::list< boost::units::angle::radian_base_unit, boost::units::list< boost::units::angle::steradian_base_unit, boost::units::dimensionless_type > > > > > > > > > >, void >").include()
#    builder.variable("kg").include()
    
    builder.global_ns.namespace('std').exclude()
    builder.build_code_creator(module_name="_chaste_project_Angiogenesis_" + module_name+"_auto")
    builder.code_creator.user_defined_directories.append( os.path.abspath('.') )
    builder.write_module(work_dir + "/dynamic/" + module_name + "_auto.cpp")

if __name__=="__main__":
    generate_wrapper_files(sys.argv)