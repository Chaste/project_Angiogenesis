/*
 * PoiseuilleImpedanceCalculator.hpp
 *
 *  Created on: 26 Feb 2015
 *      Author: chaste
 */

#ifndef PoiseuilleImpedanceCalculator_HPP_
#define PoiseuilleImpedanceCalculator_HPP_

#include "CaVascularNetwork.hpp"
#include "boost/shared_ptr.hpp"

template<unsigned DIM>
class PoiseuilleImpedanceCalculator
{

public:

	/**
	 * Constructor.
	 */
	PoiseuilleImpedanceCalculator();

	/**
	 * Destructor.
	 */
	virtual ~PoiseuilleImpedanceCalculator();

	/**
	 * Calculate impedance, Z, of all vessel segments and vessels in network using Poiseuille flow
	 * approximation:
	 *
	 * 			Z = \frac{8 \mu L}{\pi R^4},
	 *
	 * 	where \mu is viscosity, L is length and R is radius. Length is calculated within this method.
	 * 	VascularData entries "Radius" and "Viscosity" must be previously set on each segment before this
	 * 	calculation can be implemented.
	 */
	void Calculate(boost::shared_ptr<CaVascularNetwork<DIM> > vascularNetwork);

};

#endif /* PoiseuilleImpedanceCalculator_HPP_ */
