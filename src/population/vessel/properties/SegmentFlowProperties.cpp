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

template<unsigned DIM>
SegmentFlowProperties<DIM>::SegmentFlowProperties() : AbstractVesselNetworkComponentFlowProperties<DIM>(),
    mHaematocrit(0.0),
    mFlowRate(0.0*unit::unit_flow_rate),
    mImpedance(0.0*unit::unit_flow_impedance),
    mViscosity(0.0*unit::poiseuille),
    mWallShearStress(0.0*unit::pascals),
    mStimulus(0.0*unit::reciprocal_seconds)
{
}

template<unsigned DIM>
SegmentFlowProperties<DIM>::~SegmentFlowProperties()
{
}

template<unsigned DIM>
double SegmentFlowProperties<DIM>::GetHaematocrit() const
{
    return mHaematocrit;
}

template<unsigned DIM>
units::quantity<unit::dimensionless> SegmentFlowProperties<DIM>::GetDimensionalHaematocrit() const
{
    return mHaematocrit;
}

template<unsigned DIM>
double SegmentFlowProperties<DIM>::GetHaematocritSI() const
{
    return mHaematocrit;
}

template<unsigned DIM>
double SegmentFlowProperties<DIM>::GetFlowRate() const
{
    return mFlowRate/(units::pow<3>(this->mReferenceLength)/this->mReferenceTime);
}

template<unsigned DIM>
units::quantity<unit::flow_rate> SegmentFlowProperties<DIM>::GetDimensionalFlowRate() const
{
    return mFlowRate;
}

template<unsigned DIM>
double SegmentFlowProperties<DIM>::GetFlowRateSI() const
{
    return mFlowRate/unit::unit_flow_rate;
}

template<unsigned DIM>
double SegmentFlowProperties<DIM>::GetImpedance() const
{
    return mImpedance/(this->mReferenceMass/(units::pow<4>(this->mReferenceLength) * this->mReferenceTime));
}

template<unsigned DIM>
units::quantity<unit::flow_impedance> SegmentFlowProperties<DIM>::GetDimensionalImpedance() const
{
    return mImpedance;
}

template<unsigned DIM>
double SegmentFlowProperties<DIM>::GetImpedanceSI() const
{
    return mImpedance/unit::unit_flow_impedance;
}

template<unsigned DIM>
double SegmentFlowProperties<DIM>::GetViscosity() const
{
    return mViscosity/(this->mReferenceMass/(this->mReferenceLength * this->mReferenceTime));
}

template<unsigned DIM>
units::quantity<unit::dynamic_viscosity> SegmentFlowProperties<DIM>::GetDimensionalViscosity() const
{
    return mViscosity;
}

template<unsigned DIM>
double SegmentFlowProperties<DIM>::GetViscositySI() const
{
    return mViscosity/unit::poiseuille;
}

template<unsigned DIM>
double SegmentFlowProperties<DIM>::GetWallShearStress() const
{
    return mWallShearStress / (this->mReferenceMass/(this->mReferenceLength*this->mReferenceTime*this->mReferenceTime));
}

template<unsigned DIM>
units::quantity<unit::pressure> SegmentFlowProperties<DIM>::GetDimensionalWallShearStress() const
{
    return mWallShearStress;
}

template<unsigned DIM>
double SegmentFlowProperties<DIM>::GetWallShearStressSI() const
{
    return mWallShearStress / unit::pascals;
}

template<unsigned DIM>
double SegmentFlowProperties<DIM>::GetGrowthStimulus() const
{
    return mStimulus/(1.0/this->mReferenceTime);
}

template<unsigned DIM>
units::quantity<unit::rate> SegmentFlowProperties<DIM>::GetDimensionalGrowthStimulus() const
{
    return mStimulus;
}

template<unsigned DIM>
double SegmentFlowProperties<DIM>::GetGrowthStimulusSI() const
{
    return mStimulus / unit::reciprocal_seconds;
}

template<unsigned DIM>
std::map<std::string, double> SegmentFlowProperties<DIM>::GetOutputData() const
{
    std::map<std::string, double> output_data;
    output_data["Segment Haematocrit"] = this->GetHaematocrit();
    output_data["Segment Flow Rate m^3/s"] = this->GetFlowRate();
    output_data["Segment Impedance kg/m^4/s"] = this->GetImpedance();
    output_data["Segment Viscosity kg/m/s"] = this->GetViscosity();
    output_data["Segment Wall Shear Stress Pa"] = this->GetWallShearStress();
    output_data["Segment Growth Stimulus s^-1"] = this->GetGrowthStimulus();
    return output_data;
}

template<unsigned DIM>
void SegmentFlowProperties<DIM>::SetHaematocrit(double haematocrit)
{
    mHaematocrit = haematocrit;
}

template<unsigned DIM>
void SegmentFlowProperties<DIM>::SetDimensionalHaematocrit(units::quantity<unit::dimensionless> haematocrit)
{
    mHaematocrit = haematocrit;
}

template<unsigned DIM>
void SegmentFlowProperties<DIM>::SetHaematocritSI(double haematocrit)
{
    mHaematocrit = haematocrit;
}

template<unsigned DIM>
void SegmentFlowProperties<DIM>::SetFlowRate(double flowRate)
{
    mFlowRate = flowRate * (units::pow<3>(this->mReferenceLength)/this->mReferenceTime);
}

template<unsigned DIM>
void SegmentFlowProperties<DIM>::SetDimensionalFlowRate(units::quantity<unit::flow_rate> flowRate)
{
    mFlowRate = flowRate;
}

template<unsigned DIM>
void SegmentFlowProperties<DIM>::SetFlowRateSI(double flowRate)
{
    mFlowRate = flowRate * unit::unit_flow_rate;
}

template<unsigned DIM>
void SegmentFlowProperties<DIM>::SetImpedance(double impedance)
{
    mImpedance = impedance * (this->mReferenceMass/(units::pow<4>(this->mReferenceLength) * this->mReferenceTime));
}

template<unsigned DIM>
void SegmentFlowProperties<DIM>::SetDimensionalImpedance(units::quantity<unit::flow_impedance> impedance)
{
    mImpedance = impedance;
}

template<unsigned DIM>
void SegmentFlowProperties<DIM>::SetImpedanceSI(double impedance)
{
    mImpedance = impedance*unit::unit_flow_impedance;
}

template<unsigned DIM>
void SegmentFlowProperties<DIM>::SetViscosity(double viscosity)
{
    mViscosity = viscosity * (this->mReferenceMass/(this->mReferenceLength * this->mReferenceTime));
}

template<unsigned DIM>
void SegmentFlowProperties<DIM>::SetDimensionalViscosity(units::quantity<unit::dynamic_viscosity> viscosity)
{
    mViscosity = viscosity;
}

template<unsigned DIM>
void SegmentFlowProperties<DIM>::SetViscositySI(double viscosity)
{
    mViscosity = viscosity * unit::poiseuille;
}

template<unsigned DIM>
void SegmentFlowProperties<DIM>::SetWallShearStress(double value)
{
    mWallShearStress = value * (this->mReferenceMass/(this->mReferenceLength*this->mReferenceTime*this->mReferenceTime)) ;
}

template<unsigned DIM>
void SegmentFlowProperties<DIM>::SetDimensionalWallShearStress(units::quantity<unit::pressure> value)
{
    mWallShearStress = value;
}

template<unsigned DIM>
void SegmentFlowProperties<DIM>::SetWallShearStressSI(double value)
{
    mWallShearStress = value * unit::pascals;
}

template<unsigned DIM>
void SegmentFlowProperties<DIM>::SetGrowthStimulus(double value)
{
    mStimulus = value*(1.0/this->mReferenceTime);
}

template<unsigned DIM>
void SegmentFlowProperties<DIM>::SetDimensionalGrowthStimulus(units::quantity<unit::rate> value)
{
    mStimulus = value;
}

template<unsigned DIM>
void SegmentFlowProperties<DIM>::SetGrowthStimulusSI(double value)
{
    mStimulus = value*unit::reciprocal_seconds;
}

// Explicit instantiation
template class SegmentFlowProperties<2> ;
template class SegmentFlowProperties<3> ;
