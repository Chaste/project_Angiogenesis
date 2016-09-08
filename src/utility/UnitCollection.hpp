#ifndef UnitCollections_hpp
#define UnitCollections_hpp

#include <map>
#include <string>
#include <sstream>
#include <boost/units/quantity.hpp>
#include <boost/units/derived_dimension.hpp>
#include <boost/units/make_scaled_unit.hpp>
#include <boost/units/scaled_base_unit.hpp>
#include <boost/units/scale.hpp>
#include <boost/units/static_rational.hpp>
#include <boost/units/units_fwd.hpp>
#include <boost/units/io.hpp>
#include <boost/units/reduce_unit.hpp>
#include <boost/units/pow.hpp>
#include <boost/units/cmath.hpp>

#include <boost/units/systems/si.hpp>
#include <boost/units/systems/si/base.hpp>
#include <boost/units/systems/si/io.hpp>
#include <boost/units/systems/si/dimensionless.hpp>
#include <boost/units/systems/si/time.hpp>
#include <boost/units/systems/si/length.hpp>
#include <boost/units/systems/si/force.hpp>
#include <boost/units/systems/si/pressure.hpp>
#include <boost/units/systems/si/mass.hpp>
#include <boost/units/systems/si/amount.hpp>

#include <boost/units/base_units/metric/minute.hpp>
#include <boost/units/base_units/metric/hour.hpp>
#include <boost/units/base_units/metric/day.hpp>
#include <boost/units/base_units/metric/mmHg.hpp>
#include <boost/units/base_units/metric/micron.hpp>

#include "Exception.hpp"

namespace units = boost::units;
namespace unit
{
    typedef units::si::dimensionless dimensionless;

    // Time
    typedef units::si::time time;
    BOOST_UNITS_STATIC_CONSTANT(seconds, units::si::time);
    BOOST_UNITS_STATIC_CONSTANT(minutes, units::metric::minute_base_unit::unit_type);
    BOOST_UNITS_STATIC_CONSTANT(hours, units::metric::hour_base_unit::unit_type);
    BOOST_UNITS_STATIC_CONSTANT(days, units::metric::day_base_unit::unit_type);

    // Rates
    typedef units::derived_dimension<units::time_base_dimension, -1>::type rate_dimension;
    typedef units::unit<rate_dimension, units::si::system> rate;
    BOOST_UNITS_STATIC_CONSTANT(per_second, rate);
    typedef units::make_scaled_unit<rate, units::scale<60, units::static_rational<-1> > >::type per_minute_type;
    BOOST_UNITS_STATIC_CONSTANT(per_minute, per_minute_type);
    typedef units::make_scaled_unit<rate, units::scale<3600, units::static_rational<-1> > >::type per_hour_type;
    BOOST_UNITS_STATIC_CONSTANT(per_hour, per_hour_type);

    // Length
    typedef units::si::length length;
    typedef units::si::area area;
    typedef units::si::volume volume;
    BOOST_UNITS_STATIC_CONSTANT(metres, units::si::length);
//    typedef units::make_scaled_unit<units::si::length , units::scale<10, units::static_rational<-6> > >::type micron_type;
    BOOST_UNITS_STATIC_CONSTANT(microns, units::metric::micron_base_unit::unit_type);

    // Force/Pressure/Stress
    typedef units::si::force force;
    BOOST_UNITS_STATIC_CONSTANT(newtons, units::si::force);
    typedef units::si::pressure pressure;
    BOOST_UNITS_STATIC_CONSTANT(pascals, units::si::pressure);
    BOOST_UNITS_STATIC_CONSTANT(mmHg, units::metric::mmHg_base_unit::unit_type);

    // Flow
    typedef units::si::dynamic_viscosity dynamic_viscosity;
    BOOST_UNITS_STATIC_CONSTANT(poiseuille, dynamic_viscosity);

    typedef units::si::velocity velocity;

    typedef units::derived_dimension<units::length_base_dimension, 3, units::time_base_dimension, -1>::type flow_rate_dimension;
    typedef units::unit<flow_rate_dimension, units::si::system> flow_rate;
    BOOST_UNITS_STATIC_CONSTANT(metre_cubed_per_second, flow_rate);

    typedef units::derived_dimension<units::mass_base_dimension, 1, units::length_base_dimension, -4, units::time_base_dimension, -1>::type flow_impedance_dimension;
    typedef units::unit<flow_impedance_dimension, units::si::system> flow_impedance;
    BOOST_UNITS_STATIC_CONSTANT(pascal_second_per_metre_cubed, flow_impedance);

    // Mass
    typedef units::si::mass mass;
    typedef units::si::amount amount;
    BOOST_UNITS_STATIC_CONSTANT(kg, units::si::mass);
    BOOST_UNITS_STATIC_CONSTANT(moles, units::si::amount);

    typedef units::derived_dimension<units::length_base_dimension, 2, units::time_base_dimension, -1>::type diffusivity_dimension;
    typedef units::unit<diffusivity_dimension, units::si::system> diffusivity;
    BOOST_UNITS_STATIC_CONSTANT(metre_squared_per_second, diffusivity);

    typedef units::derived_dimension<units::mass_base_dimension, 1, units::time_base_dimension, -1>::type mass_flow_rate_dimension;
    typedef units::unit<mass_flow_rate_dimension, units::si::system> mass_flow_rate;
    BOOST_UNITS_STATIC_CONSTANT(kg_per_second, mass_flow_rate);

    typedef units::derived_dimension<units::amount_base_dimension, 1, units::time_base_dimension, -1>::type molar_flow_rate_dimension;
    typedef units::unit<molar_flow_rate_dimension, units::si::system> molar_flow_rate;
    BOOST_UNITS_STATIC_CONSTANT(mole_per_second, molar_flow_rate);

    typedef units::derived_dimension<units::amount_base_dimension, 1, units::length_base_dimension, -3, units::time_base_dimension, -1>::type concentration_flow_rate_dimension;
    typedef units::unit<concentration_flow_rate_dimension, units::si::system> concentration_flow_rate;
    BOOST_UNITS_STATIC_CONSTANT(mole_per_metre_cubed_per_second, concentration_flow_rate);

    typedef units::derived_dimension<units::mass_base_dimension, 1, units::length_base_dimension, -2, units::time_base_dimension, -1>::type mass_flux_dimension;
    typedef units::unit<mass_flux_dimension, units::si::system> mass_flux;
    BOOST_UNITS_STATIC_CONSTANT(kg_per_metre_squared_per_second, mass_flux);

    typedef units::derived_dimension<units::amount_base_dimension, 1, units::length_base_dimension, -2, units::time_base_dimension, -1>::type molar_flux_dimension;
    typedef units::unit<molar_flux_dimension, units::si::system> molar_flux;
    BOOST_UNITS_STATIC_CONSTANT(mole_per_metre_squared_per_second, molar_flux);

    typedef units::derived_dimension<units::amount_base_dimension, 1, units::length_base_dimension, -5, units::time_base_dimension, -1>::type concentration_flux_dimension;
    typedef units::unit<concentration_flux_dimension, units::si::system> concentration_flux;
    BOOST_UNITS_STATIC_CONSTANT(mole_per_metre_pow5_per_second, concentration_flux);

    typedef units::derived_dimension<units::amount_base_dimension, 1, units::length_base_dimension, -3>::type concentration_dimension;
    typedef units::unit<concentration_dimension, units::si::system> concentration;
    BOOST_UNITS_STATIC_CONSTANT(mole_per_metre_cubed, concentration);

    typedef units::derived_dimension<units::amount_base_dimension, 1, units::mass_base_dimension, -1>::type molar_mass_dimension;
    typedef units::unit<molar_mass_dimension, units::si::system> molar_mass;
    BOOST_UNITS_STATIC_CONSTANT(mole_per_kg, molar_mass);

    typedef units::derived_dimension<units::length_base_dimension, -3>::type number_density_dimension;
    typedef units::unit<number_density_dimension, units::si::system> number_density;
    BOOST_UNITS_STATIC_CONSTANT(per_metre_cubed, number_density);

    typedef units::derived_dimension<units::amount_base_dimension, 1, units::mass_base_dimension, -1, units::length_base_dimension, -2, units::time_base_dimension, 2>::type solubility_dimension;
    typedef units::unit<solubility_dimension, units::si::system> solubility;
    BOOST_UNITS_STATIC_CONSTANT(mole_per_metre_cubed_per_pascal, solubility);
};

class UnitTester
{
    units::quantity<unit::mass> mMyMass;

public:
    UnitTester();

    ~UnitTester();


    void SetMass(units::quantity<unit::mass> inputMass);

    units::quantity<unit::mass> GetMass();
};

#endif
