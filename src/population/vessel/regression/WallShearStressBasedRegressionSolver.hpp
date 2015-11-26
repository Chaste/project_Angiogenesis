/*
 * RegressionSolver.hpp
 *
 *  Created on: Nov 26, 2015
 *      Author: anthonyconnor
 */

#ifndef WALLSHEARSTRESSBASEDREGRESSIONSOLVER_HPP_
#define WALLSHEARSTRESSBASEDREGRESSIONSOLVER_HPP_

#include "RegressionSolver.hpp"

/**
 * This class is for simulating modifications to the vessel network due to regression.
 */
template<unsigned DIM>
class WallShearStressBasedRegressionSolver
{

    /**
     * Threshold wall shear stress level, below which vessels will be removed.
     *
     * This threshold should be prescribed in units of pascals.
     */
    double mthresholdWSSLevel;

    /**
     *  Maximum time that a vessel may exist with low wall shear stress.
     *  After this amount of time a vessel is removed completely from
     *  the vessel network.
     *
     *  This time should be prescribed in units of hours.
     */
    double mmaxTimeWithLowWallShearStress;

public:

    /**
     * Constructor.
     */
    WallShearStressBasedRegressionSolver();

    /**
     * Destructor.
     */
    virtual ~WallShearStressBasedRegressionSolver();

    /**
     *  Setter for mmaxTimeWithLowWallShearStress parameter.
     */
    void SetMaximumTimeWithLowWallShearStress(double time);

    /**
     *  Setter for mthresholdWSSLevel parameter.
     */
    void SetLowWallShearStressThreshold(double threshold);

    /**
     * Increment one step in time
     */
    virtual void Increment();

};

#endif /* WALLSHEARSTRESSBASEDREGRESSIONSOLVER_HPP_ */
