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

#include "ConstantHaematocritSolver.hpp"
#include "StructuralAdaptationSolver.hpp"
#include "UnitCollections.hpp"

template<unsigned DIM>
StructuralAdaptationSolver<DIM>::StructuralAdaptationSolver() :
        AbstractStructuralAdaptationSolver<DIM>(),
        mpHaematocritCalculator(new ConstantHaematocritSolver<DIM>),
        mpFlowSolver(new FlowSolver<DIM> ),
        mpRadiusCalculator(new RadiusCalculator<DIM> ),
        mpMetabolicStimulusCalculator(new MetabolicStimulusCalculator<DIM> ),
        mpImpedanceCalculator(new VesselImpedanceCalculator<DIM> ),
        mpMechanicalStimulusCalculator(new MechanicalStimulusCalculator<DIM> ),
        mpViscosityCalculator(new ViscosityCalculator<DIM> ),
        mpWallShearStressCalculator(new WallShearStressCalculator<DIM> )
{

}

template<unsigned DIM>
StructuralAdaptationSolver<DIM>::~StructuralAdaptationSolver()
{

}

template<unsigned DIM>
boost::shared_ptr<StructuralAdaptationSolver<DIM> > StructuralAdaptationSolver<DIM>::Create()
{
    MAKE_PTR(StructuralAdaptationSolver<DIM>, pSelf);
    return pSelf;
}

template<unsigned DIM>
void StructuralAdaptationSolver<DIM>::SetRadiusCalculator(boost::shared_ptr<RadiusCalculator<DIM> > pCalculator)
{
    mpRadiusCalculator = pCalculator;
}

template<unsigned DIM>
void StructuralAdaptationSolver<DIM>::SetMetabolicStimulusCalculator(boost::shared_ptr<MetabolicStimulusCalculator<DIM> > pCalculator)
{
    mpMetabolicStimulusCalculator = pCalculator;
}

template<unsigned DIM>
void StructuralAdaptationSolver<DIM>::SetImpedanceCalculator( boost::shared_ptr<VesselImpedanceCalculator<DIM> > pCalculator)
{
    mpImpedanceCalculator = pCalculator;
}

template<unsigned DIM>
void StructuralAdaptationSolver<DIM>::SetMechanicalStimulusCalculator(boost::shared_ptr<MechanicalStimulusCalculator<DIM> > pCalculator)
{
    mpMechanicalStimulusCalculator = pCalculator;
}

template<unsigned DIM>
void StructuralAdaptationSolver<DIM>::SetViscosityCalculator(boost::shared_ptr<ViscosityCalculator<DIM> > pCalculator)
{
    mpViscosityCalculator = pCalculator;
}

template<unsigned DIM>
void StructuralAdaptationSolver<DIM>::SetWallShearStressCalculator(boost::shared_ptr<WallShearStressCalculator<DIM> > pCalculator)
{
    mpWallShearStressCalculator = pCalculator;
}

template<unsigned DIM>
void StructuralAdaptationSolver<DIM>::SetFlowSolver(boost::shared_ptr<FlowSolver<DIM> > pCalculator)
{
    mpFlowSolver = pCalculator;
}

template<unsigned DIM>
void StructuralAdaptationSolver<DIM>::SetHaematocritCalculator(boost::shared_ptr<AbstractHaematocritSolver<DIM> > pCalculator)
{
    mpHaematocritCalculator = pCalculator;
}

template<unsigned DIM>
void StructuralAdaptationSolver<DIM>::Iterate()
{
    if(!this->mpVesselNetwork)
    {
        EXCEPTION("A vessel network is required before the SA solver can be used.");
    }

    if(SimulationTime::Instance()->GetTimeStepsElapsed()==0)
    {
        UpdateFlowSolver();
    }

    mpRadiusCalculator->SetTimestep(this->GetTimeIncrement());
    mpViscosityCalculator->Calculate(this->mpVesselNetwork);
    mpImpedanceCalculator->Calculate();
    mpFlowSolver->Update(false);
    mpFlowSolver->Solve();
    mpHaematocritCalculator->Calculate(this->mpVesselNetwork);
    mpWallShearStressCalculator->Calculate(this->mpVesselNetwork);
    mpMetabolicStimulusCalculator->Calculate();
    mpMechanicalStimulusCalculator->Calculate();
    mpRadiusCalculator->Calculate(this->mpVesselNetwork);
}

template<unsigned DIM>
void StructuralAdaptationSolver<DIM>::UpdateFlowSolver()
{
    if(mpFlowSolver)
    {
        mpFlowSolver->SetVesselNetwork(this->mpVesselNetwork);
        mpFlowSolver->SetUp();
    }
}

// Explicit instantiation
template class StructuralAdaptationSolver<2> ;
template class StructuralAdaptationSolver<3> ;
