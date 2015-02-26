/*
 * SimpleFlowSolver.cpp
 *
 *  Created on: 26 Feb 2015
 *      Author: chaste
 */

#include "SimpleFlowSolver.hpp"

template<unsigned DIM>
SimpleFlowSolver<DIM>::SimpleFlowSolver()
{


}

template<unsigned DIM>
SimpleFlowSolver<DIM>::~SimpleFlowSolver()
{


}

template<unsigned DIM>
void SimpleFlowSolver<DIM>::Implement(boost::shared_ptr<CaVascularNetwork<DIM> > vascularNetwork)
{

}

// Explicit instantiation
template class SimpleFlowSolver<2>;
template class SimpleFlowSolver<3>;

