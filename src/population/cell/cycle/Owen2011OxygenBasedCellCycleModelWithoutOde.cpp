/*

 Copyright (c) 2005-2012, University of Oxford.
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

#include "Owen2011OxygenBasedCellCycleModelWithoutOde.hpp"

#include "../mutations/CancerCellMutationState.hpp"
#include "../mutations/MacrophageMutationState.hpp"
#include "../mutations/QuiescentCancerCellMutationState.hpp"
#include "ApoptoticCellProperty.hpp"
#include "CellPropertyRegistry.hpp"
#include "CellLabel.hpp"
#include "Exception.hpp"
#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "RandomNumberGenerator.hpp"
#include "MathsCustomFunctions.hpp"

Owen2011OxygenBasedCellCycleModelWithoutOde::Owen2011OxygenBasedCellCycleModelWithoutOde() :
        mMaxRandInitialPhase(0.0), mMaxPhi(1), mCurrentQuiescentDuration(0.0), mCurrentQuiescenceOnsetTime(0.0), mEnterQuiescenceOxygenConcentration(
                8.9),  // 8.9 mmHg
        mLeaveQuiescenceOxygenConcentration(9.8),   // 9.8 mmHg
        mCriticalQuiescentDuration(4000.0)               // 4000 mins
{
    mCurrentCellCyclePhase = G_ONE_PHASE;
}

void Owen2011OxygenBasedCellCycleModelWithoutOde::ResetForDivision()
{
    assert(mReadyToDivide);
    assert(mPhi == 0);
    mCurrentCellCyclePhase = G_ONE_PHASE;
    mReadyToDivide = false;
}

double Owen2011OxygenBasedCellCycleModelWithoutOde::GetCurrentQuiescentDuration()
{
    return mCurrentQuiescentDuration;
}

double Owen2011OxygenBasedCellCycleModelWithoutOde::GetCurrentQuiescenceOnsetTime()
{
    return mCurrentQuiescenceOnsetTime;
}

void Owen2011OxygenBasedCellCycleModelWithoutOde::UpdateCellCyclePhase()
{
    // mG1Duration is set when value of phi equals 1

    if (mpCell->GetMutationState()->IsType<CancerCellMutationState>())
    {
        CheckAndLabelCell();
    }

    if (mpCell->GetMutationState()->IsType<QuiescentCancerCellMutationState>())
    {
        UpdateQuiescentDuration();
    }

    if (mpCell->GetMutationState()->IsType<MacrophageMutationState>())
    {
        mCurrentCellCyclePhase = G_ZERO_PHASE;
    }

    if (mCurrentCellCyclePhase != G_ZERO_PHASE)
    {
        double oxygen = GetOxygenConcentration();
        double dt = SimulationTime::Instance()->GetTimeStep();
        mPhi = UpdatePhi(oxygen, dt);
    }
}

bool Owen2011OxygenBasedCellCycleModelWithoutOde::ReadyToDivide()
{
    UpdateCellCyclePhase();

    // Ready to divide when mPhi exceeds 1
    if (GetPhi() > mMaxPhi)
    {
        SetPhi(0.0);
        mReadyToDivide = true;
    }
    else
    {
        mReadyToDivide = false;
    }
    return mReadyToDivide;
}

double Owen2011OxygenBasedCellCycleModelWithoutOde::GetOxygenConcentration()
{
    // Get cell's oxygen concentration
    double oxygen_concentration;

    oxygen_concentration = mpCell->GetCellData()->GetItem("oxygen");

    return (oxygen_concentration);
}

AbstractCellCycleModel* Owen2011OxygenBasedCellCycleModelWithoutOde::CreateCellCycleModel()
{
    // Create a new cell-cycle model
    Owen2011OxygenBasedCellCycleModelWithoutOde* p_model = new Owen2011OxygenBasedCellCycleModelWithoutOde();

    /*
     * Set each member variable of the new cell-cycle model that inherits
     * its value from the parent.
     *
     * Note 1: some of the new cell-cycle model's member variables (namely
     * mBirthTime, mCurrentCellCyclePhase, mReadyToDivide, mDt, mpOdeSolver)
     * will already have been correctly initialized in its constructor.
     *
     * Note 2: one or more of the new cell-cycle model's member variables
     * may be set/overwritten as soon as InitialiseDaughterCell() is called on
     * the new cell-cycle model.
     *
     * Note 3: the member variable mDimension remains unset, since this cell-cycle
     * model does not need to know the spatial dimension, so if we were to call
     * SetDimension() on the new cell-cycle model an exception would be triggered;
     * hence we do not set this member variable.
     */
    p_model->SetBirthTime(mBirthTime);
    p_model->SetMinimumGapDuration(mMinimumGapDuration);
    p_model->SetStemCellG1Duration(mStemCellG1Duration);
    p_model->SetTransitCellG1Duration(mTransitCellG1Duration);
    p_model->SetSDuration(mSDuration);
    p_model->SetG2Duration(mG2Duration);
    p_model->SetMDuration(mMDuration);
    p_model->SetEnterQuiescenceOxygenConcentration(mEnterQuiescenceOxygenConcentration);
    p_model->SetLeaveQuiescenceOxygenConcentration(mLeaveQuiescenceOxygenConcentration);
    p_model->SetCriticalQuiescentDuration(mCriticalQuiescentDuration);
    p_model->SetCurrentQuiescenceOnsetTime(mCurrentQuiescenceOnsetTime);

    return p_model;
}

void Owen2011OxygenBasedCellCycleModelWithoutOde::Initialise()
{
    mpCell->SetBirthTime(SimulationTime::Instance()->GetTime());
    SetPhi(mMaxRandInitialPhase * RandomNumberGenerator::Instance()->ranf());
}

void Owen2011OxygenBasedCellCycleModelWithoutOde::InitialiseDaughterCell()
{
    mCurrentCellCyclePhase = G_ONE_PHASE;
    SetPhi(0.0);
}

void Owen2011OxygenBasedCellCycleModelWithoutOde::UpdateQuiescentDuration()
{
    assert(mpCell->GetMutationState()->IsType<QuiescentCancerCellMutationState>());
    assert(!(mpCell->HasCellProperty<ApoptoticCellProperty>()));
    assert(!mpCell->HasApoptosisBegun());

    // Get cell's oxygen concentration
    double oxygen_concentration = GetOxygenConcentration();

    if (oxygen_concentration < mLeaveQuiescenceOxygenConcentration)
    {
        // Update the duration of the current period of hypoxia
        mCurrentQuiescentDuration = (SimulationTime::Instance()->GetTime() - mCurrentQuiescenceOnsetTime) * 60;

        if (mCurrentQuiescentDuration >= mCriticalQuiescentDuration)
        {
            /*
             * This method is usually called within a CellBasedSimulation, after the CellPopulation
             * has called CellPropertyRegistry::TakeOwnership(). This means that were we to call
             * CellPropertyRegistry::Instance() here when adding the CellLabel, we would
             * be creating a new CellPropertyRegistry. In this case the CellLabel cell
             * count would be incorrect. We must therefore access the CellLabel via the
             * cell's CellPropertyCollection.
             */
            boost::shared_ptr<AbstractCellProperty> p_cell_label =
                    mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<CellLabel>();
            mpCell->AddCellProperty(p_cell_label);
        }
    }
    else
    {
        // Reset the cell's quiescent duration.
        mCurrentQuiescentDuration = 0.0;
        mCurrentQuiescenceOnsetTime = 0.0;
        mCurrentCellCyclePhase = G_ONE_PHASE;
        /*
         * This method is usually called within a CellBasedSimulation, after the CellPopulation
         * has called CellPropertyRegistry::TakeOwnership(). This means that were we to call
         * CellPropertyRegistry::Instance() here when adding the CancerCellMutationState, we would
         * be creating a new CellPropertyRegistry. In this case the CancerCellMutationState cell
         * count would be incorrect. We must therefore access the CancerCellMutationState via the
         * cell's CellPropertyCollection.
         */
        boost::shared_ptr<AbstractCellProperty> p_cell_is_cancer =
                mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<CancerCellMutationState>();
        mpCell->SetMutationState(p_cell_is_cancer);
    }
}

void Owen2011OxygenBasedCellCycleModelWithoutOde::CheckAndLabelCell()
{
    assert(mpCell->GetMutationState()->IsType<CancerCellMutationState>());
    // Get cell's oxygen concentration
    double oxygen_concentration = GetOxygenConcentration();

    if (oxygen_concentration < mEnterQuiescenceOxygenConcentration)
    {
        /*
         * This method is usually called within a CellBasedSimulation, after the CellPopulation
         * has called CellPropertyRegistry::TakeOwnership(). This means that were we to call
         * CellPropertyRegistry::Instance() here when adding the QuiescentCancerCellMutationState, we would
         * be creating a new CellPropertyRegistry. In this case the QuiescentCancerCellMutationState cell
         * count would be incorrect. We must therefore access the QuiescentCancerCellMutationState via the
         * cell's CellPropertyCollection.
         */
        boost::shared_ptr<AbstractCellProperty> p_cell_is_quiescent =
                mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<QuiescentCancerCellMutationState>();
        mpCell->SetMutationState(p_cell_is_quiescent);

        assert(mpCell->GetMutationState()->IsType<QuiescentCancerCellMutationState>());
        assert(!(mpCell->GetMutationState()->IsType<CancerCellMutationState>()));
        mCurrentQuiescenceOnsetTime = SimulationTime::Instance()->GetTime();
        mCurrentCellCyclePhase = G_ZERO_PHASE;
    }
}

double Owen2011OxygenBasedCellCycleModelWithoutOde::GetEnterQuiescenceOxygenConcentration()
{
    return mEnterQuiescenceOxygenConcentration;
}

void Owen2011OxygenBasedCellCycleModelWithoutOde::SetEnterQuiescenceOxygenConcentration(
        double enterQuiescenceOxygenConcentration)
{
    assert(enterQuiescenceOxygenConcentration >= 0.0);
    mEnterQuiescenceOxygenConcentration = enterQuiescenceOxygenConcentration;
}

double Owen2011OxygenBasedCellCycleModelWithoutOde::GetLeaveQuiescenceOxygenConcentration()
{
    return mLeaveQuiescenceOxygenConcentration;
}

void Owen2011OxygenBasedCellCycleModelWithoutOde::SetLeaveQuiescenceOxygenConcentration(
        double leaveQuiescenceOxygenConcentration)
{
    assert(leaveQuiescenceOxygenConcentration >= 0.0);
    mLeaveQuiescenceOxygenConcentration = leaveQuiescenceOxygenConcentration;
}

double Owen2011OxygenBasedCellCycleModelWithoutOde::GetCriticalQuiescentDuration()
{
    return mCriticalQuiescentDuration;
}

void Owen2011OxygenBasedCellCycleModelWithoutOde::SetCriticalQuiescentDuration(double criticalQuiescentDuration)
{
    assert(criticalQuiescentDuration >= 0.0);
    mCriticalQuiescentDuration = criticalQuiescentDuration;
}

void Owen2011OxygenBasedCellCycleModelWithoutOde::SetCurrentQuiescenceOnsetTime(double currentQuiescenceOnsetTime)
{
    assert(currentQuiescenceOnsetTime >= 0.0);
    mCurrentQuiescenceOnsetTime = currentQuiescenceOnsetTime;
}

void Owen2011OxygenBasedCellCycleModelWithoutOde::SetPhi(double phi)
{
    mPhi = phi;
}

double Owen2011OxygenBasedCellCycleModelWithoutOde::GetPhi()
{
    return mPhi;
}

double Owen2011OxygenBasedCellCycleModelWithoutOde::UpdatePhi(double oxygen, double dt)
{
    assert(!mpCell->GetMutationState()->IsType<QuiescentCancerCellMutationState>()
            && !mpCell->GetMutationState()->IsType<MacrophageMutationState>());

    // Parameter values are taken from the Owen et al. (2011) paper. Oxygen concentration is measured in mmHg.
    //It is assumed that C=0.045 corresponds to typical dimensional perivascular oxygen tension of 20mmHg.

    if (mpCell->GetMutationState()->IsType<CancerCellMutationState>())
    {
        mC0 = 1.4; // mmHg for half maximal cell cycle progression rate
        mTmin = 1600.0; // minimum cell cycle time (corresponds to fastest possible cycle rate)
    }
    else
    {
        mC0 = 3.0; // mmHg for half maximal cell cycle progression rate
        mTmin = 3000.0; // minimum cell cycle time (corresponds to fastest possible cycle rate)
    }

    mPhi = mPhi + 60.0 * dt * oxygen / (mTmin * (mC0 + oxygen));
    SetPhi(mPhi);

    return mPhi;
}

void Owen2011OxygenBasedCellCycleModelWithoutOde::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<EnterQuiescenceOxygenConcentration>" << mEnterQuiescenceOxygenConcentration
            << "</EnterQuiescenceOxygenConcentration>\n";
    *rParamsFile << "\t\t\t<LeaveQuiescenceOxygenConcentration>" << mLeaveQuiescenceOxygenConcentration
            << "</LeaveQuiescenceOxygenConcentration>\n";
    *rParamsFile << "\t\t\t<CriticalQuiescentDuration>" << mCriticalQuiescentDuration
            << "</CriticalQuiescentDuration>\n";

    // Call method on direct parent class
    AbstractSimpleCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(Owen2011OxygenBasedCellCycleModelWithoutOde)
