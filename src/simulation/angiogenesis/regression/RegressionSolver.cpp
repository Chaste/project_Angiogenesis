/*
 * RegressionSolver.cpp
 *
 *  Created on: Nov 26, 2015
 *      Author: anthonyconnor
 */

#include "RegressionSolver.hpp"

template<unsigned DIM>
RegressionSolver<DIM>::RegressionSolver() :
    mpNetwork()
{

}

template<unsigned DIM>
RegressionSolver<DIM>::~RegressionSolver()
{

}

template<unsigned DIM>
void RegressionSolver<DIM>::SetVesselNetwork(boost::shared_ptr<VascularNetwork<DIM> > pNetwork)
{
    mpNetwork = pNetwork;
}

template<unsigned DIM>
void RegressionSolver<DIM>::Increment()
{

}

// Explicit instantiation
template class RegressionSolver<2>;
template class RegressionSolver<3>;
