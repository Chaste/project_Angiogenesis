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

#ifndef _SimpleStructuralAdaptationSolver_hpp
#define _SimpleStructuralAdaptationSolver_hpp

#include "StructuralAdaptationSolver.hpp"
#include "Alarcon03HaematocritCalculator.hpp"
#include "SimpleFlowSolver.hpp"
#include "Alarcon03MetabolicStimulusCalculator.hpp"
#include "Alarcon03RadiusCalculator.hpp"
#include "Alarcon03FlowRateCalculator.hpp"
#include "Alarcon03FlowVelocityCalculator.hpp"
#include "PoiseuilleImpedanceCalculator.hpp"
#include "Alarcon03MechanicalStimulusCalculator.hpp"
#include "Alarcon03ViscosityCalculator.hpp"
#include "Alarcon03WallShearStressCalculator.hpp"
#include "CaVascularNetwork.hpp"
#include "boost/shared_ptr.hpp"

template<unsigned DIM>
class SimpleStructuralAdaptationSolver : public StructuralAdaptationSolver<DIM>
{
    
private:

    boost::shared_ptr<AbstractHaematocritCalculator<DIM> > mHaematocritCalculator;
    boost::shared_ptr<SimpleFlowSolver<DIM> > mFlowSolver;
    
    boost::shared_ptr<Alarcon03RadiusCalculator<DIM> > mRadiusCalculator;
    boost::shared_ptr<Alarcon03MetabolicStimulusCalculator<DIM> > mMetabolicStimulusCalculator;
    boost::shared_ptr<Alarcon03FlowRateCalculator<DIM> > mFlowRateCalculator;
    boost::shared_ptr<Alarcon03FlowVelocityCalculator<DIM> > mFlowVelocityCalculator;
    boost::shared_ptr<PoiseuilleImpedanceCalculator<DIM> > mImpedanceCalculator;
    boost::shared_ptr<Alarcon03MechanicalStimulusCalculator<DIM> > mMechanicalStimulusCalculator;
    boost::shared_ptr<Alarcon03ViscosityCalculator<DIM> > mViscosityCalculator;
    boost::shared_ptr<Alarcon03WallShearStressCalculator<DIM> > mWallShearStressCalculator;
    boost::shared_ptr<AbstractVesselPropertyLocalCalculator<DIM> > mDownstreamConductedStimulusCalculator;
    boost::shared_ptr<AbstractVesselPropertyLocalCalculator<DIM> > mUpstreamConductedStimulusCalculator;

    
public:
    
    /**
     *  Default constructor.
     */
    SimpleStructuralAdaptationSolver();

    /**
     *  Virtual destructor.
     */
    virtual ~SimpleStructuralAdaptationSolver();
    
    // setters
    void SetRadiusCalculator(boost::shared_ptr<Alarcon03RadiusCalculator<DIM> > pCalulcator);
    
    void SetMetabolicStimulusCalculator(boost::shared_ptr<Alarcon03MetabolicStimulusCalculator<DIM> > pCalulcator);
    
    void SetFlowRateCalculator(boost::shared_ptr<Alarcon03FlowRateCalculator<DIM> > pCalulcator);
    
    void SetFlowVelocityCalculator(boost::shared_ptr<Alarcon03FlowVelocityCalculator<DIM> > pCalulcator);
    
    void SetImpedanceCalculator(boost::shared_ptr<PoiseuilleImpedanceCalculator<DIM> > pCalulcator);
    
    void SetMechanicalStimulusCalculator(boost::shared_ptr<Alarcon03MechanicalStimulusCalculator<DIM> > pCalulcator);
    
    void SetViscosityCalculator(boost::shared_ptr<Alarcon03ViscosityCalculator<DIM> > pCalulcator);
    
    void SetWallShearStressCalculator(boost::shared_ptr<Alarcon03WallShearStressCalculator<DIM> > pCalulcator);
    
    void SetDownstreamConductedStimulusCalculator(boost::shared_ptr<AbstractVesselPropertyLocalCalculator<DIM> > pCalulcator);
    
    void SetUpstreamConductedStimulusCalculator(boost::shared_ptr<AbstractVesselPropertyLocalCalculator<DIM> > pCalulcator);
    
    void SetNodePressureCalculator(boost::shared_ptr<SimpleFlowSolver<DIM> > pCalulcator);
    
    void SetHaematocritCalculator(boost::shared_ptr<AbstractHaematocritCalculator<DIM> > pCalulcator);
    
    // method for performing the Calculation
    virtual void Iterate(boost::shared_ptr<CaVascularNetwork<DIM> > vascularNetwork);
    
};

#endif
