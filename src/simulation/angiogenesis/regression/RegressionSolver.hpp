/*
 * RegressionSolver.hpp
 *
 *  Created on: Nov 26, 2015
 *      Author: anthonyconnor
 */

#ifndef REGRESSIONSOLVER_HPP_
#define REGRESSIONSOLVER_HPP_

#include "SmartPointers.hpp"
#include "CaVascularNetwork.hpp"

/**
 * This class is for simulating modifications to the vessel network due to regression.
 */
template<unsigned DIM>
class RegressionSolver
{

    /**
     * The vessel network
     */
    boost::shared_ptr<CaVascularNetwork<DIM> > mpNetwork;

public:

    /**
     * Constructor.
     */
    RegressionSolver();

    /**
     * Destructor.
     */
    virtual ~RegressionSolver();

    /**
     * Set the vessel network
     * @param pNetwork the vessel network
     */
    void SetVesselNetwork(boost::shared_ptr<CaVascularNetwork<DIM> > pNetwork);

    /**
     * Increment one step in time
     */
    virtual void Increment();

};

#endif /* REGRESSIONSOLVER_HPP_ */
