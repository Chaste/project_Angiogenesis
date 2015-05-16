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

#include "SegmentFlowProperties.hpp"

SegmentFlowProperties::SegmentFlowProperties() :
    mHaematocrit(0.45),
    mFlowRate(0.0),
    mImpedance(0.0),
    mViscosity(0.0),
    mWallShearStress(0.0),
    mMechanicalStimulus(0.0),
    mMetabolicStimulus(0.0),
    mUpstreamConductedStimulus(0.0),
    mDownstreamConductedStimulus(0.0),
    mShrinkingStimulus(0.0)
{
}

SegmentFlowProperties::~SegmentFlowProperties()
{
}

double SegmentFlowProperties::GetHaematocrit() const
{
    return mHaematocrit;
}

double SegmentFlowProperties::GetFlowRate() const
{
    return mFlowRate;
}

double SegmentFlowProperties::GetImpedance() const
{
    return mImpedance;
}

double SegmentFlowProperties::GetViscosity() const
{
    return mViscosity;
}

double SegmentFlowProperties::GetWallShearStress() const
{
    return mWallShearStress;
}

double SegmentFlowProperties::GetMechanicalStimulus() const
{
    return mMechanicalStimulus;
}

double SegmentFlowProperties::GetMetabolicStimulus() const
{
    return mMetabolicStimulus;
}

double SegmentFlowProperties::GetUpstreamConductedStimulus() const
{
    return mUpstreamConductedStimulus;
}

double SegmentFlowProperties::GetDownstreamConductedStimulus() const
{
    return mDownstreamConductedStimulus;
}

double SegmentFlowProperties::GetShrinkingStimulus() const
{
    return mShrinkingStimulus;
}

std::map<std::string, double> SegmentFlowProperties::GetVtkData() const
{
    std::map<std::string, double> vtk_data;
    vtk_data["Haematocrit"] = GetHaematocrit();
    vtk_data["Flow Rate"] = GetFlowRate();
    vtk_data["Impedance"] = GetImpedance();
    vtk_data["Viscosity"] = GetViscosity();
    vtk_data["Wall Shear Stress"] = GetWallShearStress();
    vtk_data["Mechanical Stimulus"] = GetMechanicalStimulus();
    vtk_data["Metabolic Stimulus"] = GetMetabolicStimulus();
    vtk_data["Upstream Conducted Stimulus"] = GetUpstreamConductedStimulus();
    vtk_data["Downstream Conducted Stimulus"] = GetDownstreamConductedStimulus();
    vtk_data["Shrinking Stimulus"] = GetShrinkingStimulus();
    return vtk_data;
}

void SegmentFlowProperties::SetHaematocrit(double haematocrit)
{
    mHaematocrit = haematocrit;
}

void SegmentFlowProperties::SetFlowRate(double flowRate)
{
    mFlowRate = flowRate;
}

void SegmentFlowProperties::SetImpedance(double impedance)
{
    mImpedance = impedance;
}

void SegmentFlowProperties::SetViscosity(double viscosity)
{
    mViscosity = viscosity;
}

void SegmentFlowProperties::SetWallShearStress(double value)
{
    mWallShearStress = value;
}

void SegmentFlowProperties::SetMechanicalStimulus(double value)
{
    mMechanicalStimulus = value;
}

void SegmentFlowProperties::SetMetabolicStimulus(double value)
{
    mMetabolicStimulus = value;
}

void SegmentFlowProperties::SetUpstreamConductedStimulus(double value)
{
    mUpstreamConductedStimulus = value;
}

void SegmentFlowProperties::SetDownstreamConductedStimulus(double value)
{
    mDownstreamConductedStimulus = value;
}

void SegmentFlowProperties::SetShrinkingStimulus(double value)
{
    mShrinkingStimulus = value;
}
