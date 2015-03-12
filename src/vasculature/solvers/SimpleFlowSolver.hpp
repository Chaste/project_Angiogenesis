/*
 * SimpleFlowSolver.hpp
 *
 *  Created on: 26 Feb 2015
 *      Author: chaste
 */

#ifndef SIMPLEFLOWSOLVER_HPP_
#define SIMPLEFLOWSOLVER_HPP_

#include <boost/shared_ptr.hpp>
#include "CaVascularNetwork.hpp"

template<unsigned DIM>
class SimpleFlowSolver
{

public:

	/**
	 * Constructor.
	 */
	SimpleFlowSolver();

	/**
	 * Destructor.
	 */
	~SimpleFlowSolver();

	/**
	 * Implement flow solver;
	 */
	void Implement(boost::shared_ptr<CaVascularNetwork<DIM> > vascularNetwork);

};

#endif /* SIMPLEFLOWSOLVER_HPP_ */
