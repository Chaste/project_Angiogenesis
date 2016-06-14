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
    mFlowRate(0.0*units::pow<3>(unit::metres)/unit::seconds),
    mImpedance(0.0*unit::kg/(units::pow<4>(unit::metres)*unit::seconds)),
    mViscosity(0.0*unit::kg/(unit::metres*unit::seconds)),
    mWallShearStress(0.0*unit::pascals),
    mStimulus(0.0*unit::reciprocal_seconds)
{
}

SegmentFlowProperties::~SegmentFlowProperties()
{
}

units::quantity<unit::dimensionless> SegmentFlowProperties::GetHaematocrit() const
{
    return mHaematocrit;
}

units::quantity<unit::flow_rate> SegmentFlowProperties::GetFlowRate() const
{
    return mFlowRate;
}

units::quantity<unit::flow_impedance> SegmentFlowProperties::GetImpedance() const
{
    return mImpedance;
}

units::quantity<unit::dynamic_viscosity> SegmentFlowProperties::GetViscosity() const
{
    return mViscosity;
}

units::quantity<unit::pressure> SegmentFlowProperties::GetWallShearStress() const
{
    return mWallShearStress;
}

units::quantity<unit::rate> SegmentFlowProperties::GetStimulus() const
{
    return mStimulus;
}

std::map<std::string, double> SegmentFlowProperties::GetVtkData() const
{
    std::map<std::string, double> vtk_data;
    vtk_data["Haematocrit"] = GetHaematocrit();
    vtk_data["Flow Rate m3/s"] = GetFlowRate()/(units::pow<3>(unit::metres)/unit::seconds);
    vtk_data["Impedance kg/m4/s"] = GetImpedance()/(unit::kg/(units::pow<4>(unit::metres)*unit::seconds));
    vtk_data["Viscosity kg/m/s"] = GetViscosity()/(unit::kg/(unit::metres*unit::seconds));
    vtk_data["Wall Shear Stress Pa"] = GetWallShearStress()/(unit::pascals);
    vtk_data["Growth Stimulus s-1"] = GetStimulus()/(unit::reciprocal_seconds);
    return vtk_data;
}

void SegmentFlowProperties::SetHaematocrit(units::quantity<unit::dimensionless> haematocrit)
{
    mHaematocrit = haematocrit;
}

void SegmentFlowProperties::SetFlowRate(units::quantity<unit::flow_rate> flowRate)
{
    mFlowRate = flowRate;
}

void SegmentFlowProperties::SetImpedance(units::quantity<unit::flow_impedance> impedance)
{
    mImpedance = impedance;
}

void SegmentFlowProperties::SetViscosity(units::quantity<unit::dynamic_viscosity> viscosity)
{
    mViscosity = viscosity;
}

void SegmentFlowProperties::SetWallShearStress(units::quantity<unit::pressure> value)
{
    mWallShearStress = value;
}

void SegmentFlowProperties::SetStimulus(units::quantity<unit::rate> value)
{
    mStimulus = value;
}
