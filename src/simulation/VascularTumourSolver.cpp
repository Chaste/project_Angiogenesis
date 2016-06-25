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

#include <boost/lexical_cast.hpp>
#include "UblasIncludes.hpp"
#include "VesselSegment.hpp"
#include "VesselNode.hpp"
#include "VascularTumourSolver.hpp"
#include "VtkVesselNetworkWriter.hpp"

template<unsigned DIM>
VascularTumourSolver<DIM>::VascularTumourSolver() :
        mpNetwork(),
        mOutputFrequency(1),
        mpOutputFileHandler(),
        mHybridSolvers(),
        mpStructuralAdaptationSolver(),
        mpAngiogenesisSolver(),
        mpRegressionSolver()
{

}

template<unsigned DIM>
VascularTumourSolver<DIM>::~VascularTumourSolver()
{

}

template<unsigned DIM>
boost::shared_ptr<VascularTumourSolver<DIM> > VascularTumourSolver<DIM>::Create()
{
    MAKE_PTR(VascularTumourSolver<DIM>, pSelf);
    return pSelf;
}

template<unsigned DIM>
void VascularTumourSolver<DIM>::SetVesselNetwork(boost::shared_ptr<VesselNetwork<DIM> > pNetwork)
{
    mpNetwork = pNetwork;
}

template<unsigned DIM>
void VascularTumourSolver<DIM>::AddHybridSolver(boost::shared_ptr<AbstractHybridSolver<DIM> > pHybridSolver)
{
    mHybridSolvers.push_back(pHybridSolver);
}

template<unsigned DIM>
std::vector<boost::shared_ptr<AbstractHybridSolver<DIM> > > VascularTumourSolver<DIM>::GetHybridSolvers()
{
    return mHybridSolvers;
}
template<unsigned DIM>
void VascularTumourSolver<DIM>::Increment()
{
    unsigned num_steps = SimulationTime::Instance()->GetTimeStepsElapsed();

    // If there is a structural adaptation or flow problem solve it
    if(mpStructuralAdaptationSolver)
    {
        mpStructuralAdaptationSolver->UpdateFlowSolver();
        mpStructuralAdaptationSolver->Solve();
    }

    // If there are PDEs solve them
    if(mHybridSolvers.size()>0)
    {
        for(unsigned idx=0; idx<mHybridSolvers.size(); idx++)
        {
            mHybridSolvers[idx]->Update();
            mHybridSolvers[idx]->SetFileName("/" + mHybridSolvers[idx]->GetLabel() +"_solution_" + boost::lexical_cast<std::string>(num_steps)+".vti");

            // Transfer PDE solutions for coupled problems
//            if(idx>0 and mHybridSolvers[idx]->GetPde())
//            {
//                for(unsigned jdx=0; jdx<mHybridSolvers[idx]->GetPde()->GetDiscreteSources().size(); jdx++)
//                {
//                    if(mHybridSolvers[idx]->GetPde()->GetDiscreteSources()[jdx]->GetType()==SourceType::SOLUTION)
//                    {
//                        mHybridSolvers[idx]->GetPde()->GetDiscreteSources()[jdx]->SetSolution(mHybridSolvers[idx-1]->GetSolutionAtGridPoints());
//                    }
//                }
//            }
            if(mOutputFrequency > 0 && num_steps % mOutputFrequency == 0)
            {
                mHybridSolvers[idx]->SetWriteSolution(true);
                mHybridSolvers[idx]->Solve();
            }
            else
            {
                mHybridSolvers[idx]->SetWriteSolution(false);
                mHybridSolvers[idx]->Solve();
            }
        }
    }

    // Do angiogenesis if the is a network and solver
    if(this->mpNetwork && mpAngiogenesisSolver)
    {
        mpAngiogenesisSolver->Increment();
    }

    // Do regression if the is a network and solver
    if(this->mpNetwork && mpRegressionSolver)
    {
        mpRegressionSolver->Increment();
    }

    // Manage vessel network output
    if(this->mpNetwork)
    {
        mpNetwork->UpdateAll();
        if(mOutputFrequency > 0 && num_steps % mOutputFrequency == 0)
        {
        	boost::shared_ptr<VtkVesselNetworkWriter<DIM> > p_network_writer = VtkVesselNetworkWriter<DIM>::Create();
    		p_network_writer->SetVesselNetwork(mpNetwork);
    		p_network_writer->SetFileName(mpOutputFileHandler->GetOutputDirectoryFullPath() + "/VesselNetwork_inc_" +
    				boost::lexical_cast<std::string>(num_steps+1)+".vtp");
    		p_network_writer->Write();
        }
    }
}

template<unsigned DIM>
void VascularTumourSolver<DIM>::Run()
{
    if(this->mpNetwork)
    {
    	boost::shared_ptr<VtkVesselNetworkWriter<DIM> > p_network_writer = VtkVesselNetworkWriter<DIM>::Create();
        mpNetwork->UpdateAll(true);
		p_network_writer->SetVesselNetwork(mpNetwork);
		p_network_writer->SetFileName(mpOutputFileHandler->GetOutputDirectoryFullPath() + "/VesselNetwork_inc_0.vtp");
		p_network_writer->Write();
    }

    Setup();
    while(!SimulationTime::Instance()->IsFinished())
    {
        Increment();
        SimulationTime::Instance()->IncrementTimeOneStep();
    }
}

template<unsigned DIM>
void VascularTumourSolver<DIM>::SetAngiogenesisSolver(boost::shared_ptr<AngiogenesisSolver<DIM> > pAngiogenesisSolver)
{
    mpAngiogenesisSolver = pAngiogenesisSolver;
}

template<unsigned DIM>
void VascularTumourSolver<DIM>::SetOutputFileHandler(boost::shared_ptr<OutputFileHandler> pFileHandler)
{
    mpOutputFileHandler = pFileHandler;
}

template<unsigned DIM>
void VascularTumourSolver<DIM>::SetOutputFrequency(unsigned frequency)
{
    mOutputFrequency = frequency;
}

template<unsigned DIM>
void VascularTumourSolver<DIM>::SetupFromModifier(AbstractCellPopulation<DIM,DIM>& rCellPopulation, const std::string& rDirectory)
{
    // Set up an output file handler
    mpOutputFileHandler = boost::shared_ptr<OutputFileHandler>(new OutputFileHandler(rDirectory));

    // Set up all the hybrid solvers
    for(unsigned idx=0; idx<mHybridSolvers.size(); idx++)
    {
        mHybridSolvers[idx]->SetCellPopulation(rCellPopulation);
    }

    Setup();
}

template<unsigned DIM>
void VascularTumourSolver<DIM>::Setup()
{
    // Set up all the hybrid solvers
    for(unsigned idx=0; idx<mHybridSolvers.size(); idx++)
    {
        mHybridSolvers[idx]->SetFileHandler(mpOutputFileHandler);
        if(mpNetwork)
        {
            mHybridSolvers[idx]->SetVesselNetwork(mpNetwork);
        }
        mHybridSolvers[idx]->Setup();
    }

    // Set up the flow and structural adaptation solvers
    if(mpStructuralAdaptationSolver)
    {
        mpStructuralAdaptationSolver->SetVesselNetwork(mpNetwork);
    }

    if(mpRegressionSolver)
    {
        mpRegressionSolver->SetVesselNetwork(mpNetwork);
    }
}

template<unsigned DIM>
void VascularTumourSolver<DIM>::SetStructuralAdaptationSolver(boost::shared_ptr<StructuralAdaptationSolver<DIM> > pStructuralAdaptationSolver)
{
    mpStructuralAdaptationSolver = pStructuralAdaptationSolver;
}

template<unsigned DIM>
void VascularTumourSolver<DIM>::UpdateCellData(std::vector<std::string> labels)
{
    if(labels.size()==0)
    {
        //update everything
        for(unsigned jdx=0; jdx<mHybridSolvers.size(); jdx++)
        {
            mHybridSolvers[jdx]->UpdateCellData();
        }
    }

    for(unsigned idx=0; idx<labels.size(); idx++)
    {
        for(unsigned jdx=0; jdx<mHybridSolvers.size(); jdx++)
        {
            if(labels[idx]==mHybridSolvers[idx]->GetPde()->GetVariableName())
            {
                mHybridSolvers[jdx]->UpdateCellData();
            }
        }
    }
}

template<unsigned DIM>
void VascularTumourSolver<DIM>::SetRegressionSolver(boost::shared_ptr<RegressionSolver<DIM> > pRegressionSolver)
{
    mpRegressionSolver = pRegressionSolver;
}

// Explicit instantiation
template class VascularTumourSolver<2> ;
template class VascularTumourSolver<3> ;
