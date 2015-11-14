/*

 Copyright (c) 2005-2015, University of Oxford.
 All rights reserved.

 University of Oxford means the Chancellor, Masters and Scholars of the
 University of Oxford, having an administrative office at Wellington
 Square, Oxford OX1 2JD, UK.

 This file is part of Chaste.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
 this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
 this list of conditions and the following disclaimer in the documentation
 and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
 contributors may be used to endorse or promote products derived from this
 software without specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
 OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */

#include "AbstractHybridSolver.hpp"

template<unsigned DIM>
AbstractHybridSolver<DIM>::AbstractHybridSolver()
    :   mpNetwork(),
        mpCellPopulation(NULL),
        mCellPopulationIsSet(false),
        mpPde(),
        mBoundaryConditions(),
        mpOutputFileHandler(),
        mFilename(),
        mLabel("Default"),
        mPointSolution(),
        IsSetupForSolve(false)
{

}

template<unsigned DIM>
AbstractHybridSolver<DIM>::~AbstractHybridSolver()
{

}

template<unsigned DIM>
void AbstractHybridSolver<DIM>::AddBoundaryCondition(boost::shared_ptr<HybridBoundaryCondition<DIM> > pBoundaryCondition)
{
    mBoundaryConditions.push_back(pBoundaryCondition);
}

template<unsigned DIM>
boost::shared_ptr<HybridLinearEllipticPde<DIM, DIM> > AbstractHybridSolver<DIM>::GetPde()
{
    if(!mpPde)
    {
        EXCEPTION("A pde has not been set.");
    }
    return mpPde;
}

template<unsigned DIM>
void AbstractHybridSolver<DIM>::SetCellPopulation(AbstractCellPopulation<DIM>& rCellPopulation)
{
    mpCellPopulation = &rCellPopulation;
    mCellPopulationIsSet = true;
}

template<unsigned DIM>
std::vector<double> AbstractHybridSolver<DIM>::GetPointSolution()
{
    return mPointSolution;
}

template<unsigned DIM>
void AbstractHybridSolver<DIM>::SetFileHandler(boost::shared_ptr<OutputFileHandler> pOutputFileHandler)
{
    mpOutputFileHandler = pOutputFileHandler;
}

template<unsigned DIM>
void AbstractHybridSolver<DIM>::SetFileName(const std::string& rFilename)
{
    mFilename = rFilename;
}

template<unsigned DIM>
void AbstractHybridSolver<DIM>::SetLabel(const std::string& label)
{
    mLabel = label;
}

template<unsigned DIM>
const std::string& AbstractHybridSolver<DIM>::GetLabel()
{
    return mLabel;
}

template<unsigned DIM>
void AbstractHybridSolver<DIM>::SetPde(boost::shared_ptr<HybridLinearEllipticPde<DIM, DIM> > pPde)
{
    mpPde = pPde;
}

template<unsigned DIM>
void AbstractHybridSolver<DIM>::SetVesselNetwork(boost::shared_ptr<CaVascularNetwork<DIM> > pNetwork)
{
    mpNetwork = pNetwork;
}

// Explicit instantiation
template class AbstractHybridSolver<2> ;
template class AbstractHybridSolver<3> ;
