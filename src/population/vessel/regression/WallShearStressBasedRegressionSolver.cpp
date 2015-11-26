/*
 * RegressionSolver.cpp
 *
 *  Created on: Nov 26, 2015
 *      Author: anthonyconnor
 */

#include "WallShearStressBasedRegressionSolver.hpp"
#include <vector>

template<unsigned DIM>
WallShearStressBasedRegressionSolver<DIM>::WallShearStressBasedRegressionSolver() :
    RegressionSolver(),
    mthresholdWSSLevel(9),
    mmaxTimeWithLowWallShearStress(60)
{

}

template<unsigned DIM>
WallShearStressBasedRegressionSolver<DIM>::~WallShearStressBasedRegressionSolver()
{

}

template <int DIM, class T>
void WallShearStressBasedRegressionSolver<DIM,T>::SetMaximumTimeWithLowWallShearStress(double time)
{
    assert(time > 0);
    mmaxTimeWithLowWallShearStress = time;
}

template <int DIM, class T>
void WallShearStressBasedRegressionSolver<DIM,T>::SetLowWallShearStressThreshold(double threshold)
{
    assert(threshold >= 0);
    mthresholdWSSLevel = threshold;
}

template<unsigned DIM>
void WallShearStressBasedRegressionSolver<DIM>::Increment()
{

    if(!this->mpNetwork)
    {
        EXCEPTION("The regression solver needs an initial vessel network");
    }

    std::vector<boost::shared_ptr<CaVessel<DIM> > > vessels = mpNetwork->GetVessels();
    unsigned numberOfVessels = vessels.size();

    for(std::vector<boost::shared_ptr<CaVessel<DIM> > >::iterator v_it = vessels.begin(); v_it != vessels.end(); v_it++)
    {
        assert(v_it->GetSegment(0)->GetFlowProperties()->GetWallShearStress() >= 0);

        // if wall shear stress of vessel is below threshold then start regression timer unless it has already been started
        if (v_it->GetSegment(0)->GetFlowProperties()->GetWallShearStress() < mthresholdWSSLevel)
        {

            if (!(v_it->HasRegressionTimerStarted()) && !(v_it->VesselHasRegressed()))
            {
                // increment time that the vessel has had low wall shear stress
                v_it->SetTimeUntilRegression(mmaxTimeWithLowWallShearStress);
            }

        }
        else // otherwise rescue vessel
        {
            // wall shear stress above threshold so vessel is not regressing
            v_it->ResetRegressionTimer();
        }
    }

    assert(vessels.size() == numberOfVessels);

    // iterate through all vessels and if regression flag is true then remove from the network
    for(std::vector<boost::shared_ptr<CaVessel<DIM> > >::iterator v_it = vessels.begin(); v_it != vessels.end(); v_it++)
    {
        if (v_it->VesselHasRegressed())
        {
            mpNetwork->RemoveVessel(*v_it, true); // should this flag be true ... method needs better documentation
            assert(vessels.size() == numberOfVessels); // sanity checking
        }
    }

}

// Explicit instantiation
template class WallShearStressBasedRegressionSolver<2>;
template class WallShearStressBasedRegressionSolver<3>;
