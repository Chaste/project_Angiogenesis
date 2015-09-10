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
        mpCellPopulation(),
        mpPde(),
        mBoundaryConditionType(BoundaryConditionType::LINE),
        mBoundaryConditionValue(0.0),
        mBoundaryConditionSource(BoundaryConditionSource::USER),
        mBoundaryConditionName(""),
        mpBoundaryCondition(),
        mpInterfaceCondition(),
        mWorkingDirectory(),
        mFilename()
{

}

template<unsigned DIM>
AbstractHybridSolver<DIM>::~AbstractHybridSolver()
{

}

template<unsigned DIM>
void AbstractHybridSolver<DIM>::SetBoundaryConditionType(BoundaryConditionType::Value boundaryType)
{
    mBoundaryConditionType = boundaryType;
}

template<unsigned DIM>
void AbstractHybridSolver<DIM>::SetBoundaryConditionSource(BoundaryConditionSource::Value boundarySource)
{
    mBoundaryConditionSource = boundarySource;
}

template<unsigned DIM>
void AbstractHybridSolver<DIM>::SetBoundaryConditionValue(double boundaryConditionValue)
{
    mBoundaryConditionValue = boundaryConditionValue;
}

template<unsigned DIM>
void AbstractHybridSolver<DIM>::SetBoundaryConditionName(const std::string& rBoundaryConditionName)
{
    mBoundaryConditionName = rBoundaryConditionName;
}

template<unsigned DIM>
void AbstractHybridSolver<DIM>::SetFileName(const std::string& filename)
{
    mFilename = filename;
}

template<unsigned DIM>
void AbstractHybridSolver<DIM>::SetWorkingDirectory(const std::string& rDirectory)
{
    mWorkingDirectory = rDirectory;
}

template<unsigned DIM>
vtkSmartPointer<vtkImageData> AbstractHybridSolver<DIM>::GetSolution()
{
    return vtkSmartPointer<vtkImageData>::New();
}

template<unsigned DIM>
void AbstractHybridSolver<DIM>::SetCellPopulation(boost::shared_ptr<SimpleCellPopulation<DIM> > pCellPopulation)
{
    mpCellPopulation = pCellPopulation;
}

template<unsigned DIM>
void AbstractHybridSolver<DIM>::SetDomainBoundaryCondition(boost::shared_ptr<AbstractBoundaryCondition<DIM> > pBoundaryCondition)
{
    mpBoundaryCondition = pBoundaryCondition;
}

template<unsigned DIM>
void AbstractHybridSolver<DIM>::SetInterfaceBoundaryCondition(boost::shared_ptr<AbstractBoundaryCondition<DIM> > pBoundaryCondition)
{
    mpInterfaceCondition = pBoundaryCondition;
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

template<unsigned DIM>
void AbstractHybridSolver<DIM>::Solve(bool writeSolution)
{
    if (writeSolution)
    {
        Write();
    }
}

template<unsigned DIM>
void AbstractHybridSolver<DIM>::UpdateSolution(std::map<std::string, std::vector<double> >& data)
{

}

template<unsigned DIM>
void AbstractHybridSolver<DIM>::Write()
{

}

// Explicit instantiation
template class AbstractHybridSolver<2> ;
template class AbstractHybridSolver<3> ;
