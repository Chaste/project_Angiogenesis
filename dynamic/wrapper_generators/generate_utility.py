#!/usr/bin/env python

"""
This scipt automatically generates Python bindings using a rule based approach
"""
import sys
from pyplusplus import module_builder
from pyplusplus.module_builder import call_policies
from pygccxml import parser

def update_builder(builder):
    
    # Manually build up the boost units names, start with the SI system
    prefix = "boost::units::list<boost::units::"
    solid_angle_list = prefix + "angle::steradian_base_unit, boost::units::dimensionless_type>"
    angle_list = prefix + "angle::radian_base_unit, "+ solid_angle_list +">"
    luminosity_list = prefix + "si::candela_base_unit, "+ angle_list +">"
    amount_list = prefix + "si::mole_base_unit, "+ luminosity_list +">"
    temperature_list = prefix + "si::kelvin_base_unit, "+ amount_list +">"
    current_list = prefix + "si::ampere_base_unit, "+ temperature_list +">"
    time_list = prefix + "si::second_base_unit, "+ current_list +">"
    #mass_list = prefix + "scaled_base_unit<boost::units::cgs::gram_base_unit, boost::units::scale<10, boost::units::static_rational<3> > >, "+ time_list +">"
    mass_list = prefix + "scaled_base_unit<boost::units::cgs::gram_base_unit, boost::units::scale<10, static_rational<3> > >, "+ time_list +">"
    length_list = prefix + "si::meter_base_unit, "+ mass_list +">"
    si_homogeneous_system = "boost::units::homogeneous_system<"+ length_list +">"

    # Next each unit and quantity
    # The building blocks
    mass_pow_1 = "boost::units::list< boost::units::dim< boost::units::mass_base_dimension, boost::units::static_rational<1, 1> >, boost::units::dimensionless_type>"
    length_pow_1 = "boost::units::list< boost::units::dim< boost::units::length_base_dimension, boost::units::static_rational< 1, 1 > >, boost::units::dimensionless_type >"
    time_pow_1 = "boost::units::list<boost::units::dim<boost::units::time_base_dimension, boost::units::static_rational<1, 1> >, boost::units::dimensionless_type>"
    time_pow_neg_1 = "boost::units::list< boost::units::dim< boost::units::time_base_dimension, boost::units::static_rational< -1, 1 > >, boost::units::dimensionless_type >"
    
    dimensionless_unit = "unit< boost::units::dimensionless_type, " +  si_homogeneous_system +", void >"

    include_classes = ["UnitTester", 
                       "ParameterCollection", 
#                       "BaseParameterInstance", 
#                       "DimensionalSimulationTime",
#                       "LengthParameterInstance",
#                       "MassParameterInstance",
#                       "TimeParameterInstance",
#                       "PressureParameterInstance",
                       "ViscosityParameterInstance"]
    
    for eachClass in include_classes:
        builder.class_(eachClass).include()

    #builder.class_("kg_instance_t<true>").include()
    #builder.class_("metre_cubed_per_second_instance_t< true >").include()
#    builder.class_(dimensionless_unit).include()
#    builder.class_("DimensionlessUnit").include()
    
    # can we find all the boost units classes 
    helpers = builder.classes(lambda decl: decl.name.startswith('unit'))
    helpers.include()
    
#     helpers = builder.classes(lambda decl: decl.name.startswith('quantity'))
#     helpers.include()
    
    
#     PressureUnit@::pyplusplus::aliases::PressureUnit
#     builder.print_declarations()
    
#    builder.class_("boost::units::unit< boost::units::list< boost::units::dim< boost::units::mass_base_dimension, boost::units::static_rational< 1, 1 > >, boost::units::dimensionless_type >, boost::units::homogeneous_system< boost::units::list< boost::units::si::meter_base_unit, boost::units::list< boost::units::scaled_base_unit< boost::units::cgs::gram_base_unit, boost::units::scale< 10, static_rational< 3 > > >, boost::units::list< boost::units::si::second_base_unit, boost::units::list< boost::units::si::ampere_base_unit, boost::units::list< boost::units::si::kelvin_base_unit, boost::units::list< boost::units::si::mole_base_unit, boost::units::list< boost::units::si::candela_base_unit, boost::units::list< boost::units::angle::radian_base_unit, boost::units::list< boost::units::angle::steradian_base_unit, boost::units::dimensionless_type > > > > > > > > > >, void >").include()
    builder.variable("kg").include()

#    builder.class_('::boost::units::quantity< boost::units::unit< boost::units::list< boost::units::dim< boost::units::mass_base_dimension, boost::units::static_rational< 1, 1 > >, boost::units::dimensionless_type >, boost::units::homogeneous_system< boost::units::list< boost::units::si::meter_base_unit, boost::units::list< boost::units::scaled_base_unit< boost::units::cgs::gram_base_unit, boost::units::scale< 10, boost::units::static_rational< 3 > > >, boost::units::list< boost::units::si::second_base_unit, boost::units::list< boost::units::si::ampere_base_unit, boost::units::list< boost::units::si::kelvin_base_unit, boost::units::list< boost::units::si::mole_base_unit, boost::units::list< boost::units::si::candela_base_unit, boost::units::list< boost::units::angle::radian_base_unit, boost::units::list< boost::units::angle::steradian_base_unit, boost::units::dimensionless_type > > > > > > > > > >, void >, double >').add_registration_code('def(float() * bp::self)')
#    builder.class_('::boost::units::quantity< boost::units::unit< boost::units::list< boost::units::dim< boost::units::mass_base_dimension, boost::units::static_rational< 1, 1 > >, boost::units::dimensionless_type >, boost::units::homogeneous_system< boost::units::list< boost::units::si::meter_base_unit, boost::units::list< boost::units::scaled_base_unit< boost::units::cgs::gram_base_unit, boost::units::scale< 10, boost::units::static_rational< 3 > > >, boost::units::list< boost::units::si::second_base_unit, boost::units::list< boost::units::si::ampere_base_unit, boost::units::list< boost::units::si::kelvin_base_unit, boost::units::list< boost::units::si::mole_base_unit, boost::units::list< boost::units::si::candela_base_unit, boost::units::list< boost::units::angle::radian_base_unit, boost::units::list< boost::units::angle::steradian_base_unit, boost::units::dimensionless_type > > > > > > > > > >, void >, double >').add_registration_code('def(bp::self * float())')
    
    return builder