/*
 * SimpleFlowSolver.hpp
 *
 *  Created on: 26 Feb 2015
 *      Author: chaste
 */

#ifndef SIMPLEFLOWSOLVER_HPP_
#define SIMPLEFLOWSOLVER_HPP_

#include "CaVascularNetwork.hpp"
#include "boost/shared_ptr.hpp"

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
	virtual ~SimpleFlowSolver();

	/**
	 * Implement flow solver;
	 */
	void Implement(boost::shared_ptr<CaVascularNetwork<DIM> > vascularNetwork);

};

#endif /* SIMPLEFLOWSOLVER_HPP_ */
