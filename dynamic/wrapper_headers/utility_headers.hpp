#include "BaseUnits.hpp"
#include "LengthParameterInstance.hpp"
#include "MassParameterInstance.hpp"
#include "TimeParameterInstance.hpp"
#include "PressureParameterInstance.hpp"
#include "ViscosityParameterInstance.hpp"
#include "UnitCollection.hpp"
#include "ParameterCollection.hpp"
#include "BaseParameterInstance.hpp"

// Typdef in this namespace so that pyplusplus uses the nicer typedef'd name for the class
namespace pyplusplus{
namespace aliases{
using namespace boost::units;

// Need to break the boost units library down into more simple parts, first the SI unit system
typedef list<angle::steradian_base_unit, dimensionless_type> solid_angle_list;
typedef list<angle::radian_base_unit, solid_angle_list> angle_list;
typedef list<si::candela_base_unit, angle_list > luminosity_list;
typedef list<si::mole_base_unit, luminosity_list > amount_list;
typedef list<si::kelvin_base_unit, amount_list > temperature_list;
typedef list<si::ampere_base_unit, temperature_list > current_list;
typedef list<si::second_base_unit, current_list > time_list;
typedef list<scaled_base_unit<cgs::gram_base_unit, scale<10, static_rational<3> > >, time_list > mass_list; // i.e. kilogram
typedef list<si::meter_base_unit, mass_list > length_list;
typedef homogeneous_system< length_list> si_homogeneous_system;

// Next each unit and quantity
// The building blocks
typedef list< dim< mass_base_dimension, static_rational<1, 1> >, dimensionless_type> mass_pow_1;
typedef list< dim< length_base_dimension, static_rational< 1, 1 > >, dimensionless_type > length_pow_1;
typedef list<dim<time_base_dimension, static_rational<1, 1> >, dimensionless_type> time_pow_1;
typedef list< dim< time_base_dimension, static_rational< -1, 1 > >, dimensionless_type > time_pow_neg_1;

// Dimensionless
typedef boost::units::unit< dimensionless_type, si_homogeneous_system, void > DimensionlessUnit;
typedef quantity<DimensionlessUnit, double> DimensionlessQuantity;
//template class DimensionlessUnit;
//template class DimensionlessQuantity;

// Mass
typedef boost::units::unit< mass_pow_1, si_homogeneous_system, void> MassUnit;
typedef quantity<MassUnit, double> MassQuantity;
//template class MassUnit;
//template class MassQuantity;

// Length
typedef boost::units::unit< length_pow_1, si_homogeneous_system, void > LengthUnit;
typedef quantity<LengthUnit, double > LengthQuantity;
//template class LengthUnit;
//template class LengthQuantity;

// Time
typedef boost::units::unit<time_pow_1 , si_homogeneous_system, void> TimeUnit;
typedef quantity<TimeUnit, double > TimeQuantity;
//template class TimeUnit;
//template class TimeQuantity;

// Flow Rate
typedef boost::units::unit< list< dim< length_base_dimension, static_rational< 3, 1 > >,
        time_pow_neg_1>, si_homogeneous_system, void > FlowRateUnit;
typedef quantity<FlowRateUnit, double > FlowRateQuantity;
//template class FlowRateUnit;
//template class FlowRateQuantity;

// Viscosity
typedef boost::units::unit< list< dim< length_base_dimension, static_rational< -1, 1 > >,
        list< dim< mass_base_dimension, static_rational< 1, 1 > >,
        time_pow_neg_1 > >, si_homogeneous_system, void > ViscosityUnit;
typedef quantity<ViscosityUnit, double > ViscosityQuantity;
//template class ViscosityUnit;
//template class ViscosityQuantity;

// Impedance
typedef boost::units::unit< list< dim< length_base_dimension, static_rational< -4, 1 > >,
        list< dim< mass_base_dimension, static_rational< 1, 1 > >,
        time_pow_neg_1 > >, si_homogeneous_system, void > ImpedanceUnit;
typedef quantity<ImpedanceUnit, double > ImpedanceQuantity;
//template class ImpedanceUnit;
//template class ImpedanceQuantity;

// Pressure
typedef boost::units::unit< list< dim< length_base_dimension, static_rational< -1, 1 > >,
        list< dim< mass_base_dimension, static_rational< 1, 1 > >,
        list< dim< time_base_dimension, static_rational< -2, 1 > >,
        dimensionless_type > > >, si_homogeneous_system, void > PressureUnit;
typedef quantity<PressureUnit, double > PressureQuantity;


}
}//pyplusplus::aliases

inline int Instantiation()
{
    return  sizeof(unit::flow_rate);
}
//template class PressureUnit;
//template class PressureQuantity;

//typedef unit::kg_instance_t< true > kg;

