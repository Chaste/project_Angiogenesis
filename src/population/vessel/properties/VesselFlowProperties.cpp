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

#include "VesselFlowProperties.hpp"

template<unsigned DIM>
VesselFlowProperties<DIM>::VesselFlowProperties() : AbstractVesselNetworkComponentFlowProperties<DIM>(),
    mUndergoingRegression(false),
    mRemoveViaRegression(false),
    mRegressionTime(DBL_MAX*unit::seconds)
{
}

template<unsigned DIM>
VesselFlowProperties<DIM>::~VesselFlowProperties()
{
}

template<unsigned DIM>
units::quantity<unit::dimensionless> VesselFlowProperties<DIM>::GetDimensionalHaematocrit(const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments) const
{
    units::quantity<unit::dimensionless> haematocrit = 0.0;
    for (unsigned i = 0; i < segments.size(); i++)
    {
        haematocrit += segments[i]->GetFlowProperties()->GetDimensionalHaematocrit();
    }
    return haematocrit / double(segments.size());
}

template<unsigned DIM>
double VesselFlowProperties<DIM>::GetHaematocrit(const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments) const
{
    return GetDimensionalHaematocrit(segments);
}

template<unsigned DIM>
double VesselFlowProperties<DIM>::GetHaematocritSI(const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments) const
{
    return GetHaematocritSI(segments);
}

template<unsigned DIM>
units::quantity<unit::flow_rate> VesselFlowProperties<DIM>::GetDimensionalFlowRate(const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments) const
{
    units::quantity<unit::flow_rate> value = 0.0 * unit::unit_flow_rate;
    for (unsigned i = 0; i < segments.size(); i++)
    {
        value += segments[i]->GetFlowProperties()->GetDimensionalFlowRate();
    }
    return value / double(segments.size());
}

template<unsigned DIM>
double VesselFlowProperties<DIM>::GetFlowRate(const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments) const
{
    return GetDimensionalFlowRate(segments) /(units::pow<3>(this->mReferenceLength)/this->mReferenceTime);
}

template<unsigned DIM>
double VesselFlowProperties<DIM>::GetFlowRateSI(const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments) const
{
    return GetDimensionalFlowRate(segments) / unit::unit_flow_rate;
}

template<unsigned DIM>
units::quantity<unit::flow_impedance> VesselFlowProperties<DIM>::GetDimensionalImpedance(const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments) const
{
    units::quantity<unit::flow_impedance> value = 0.0 * unit::unit_flow_impedance;
    for (unsigned i = 0; i < segments.size(); i++)
    {
        value += segments[i]->GetFlowProperties()->GetDimensionalImpedance();
    }
    return value;
}

template<unsigned DIM>
double VesselFlowProperties<DIM>::GetImpedance(const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments) const
{
    return GetDimensionalImpedance(segments) /(this->mReferenceMass/(units::pow<4>(this->mReferenceLength) * this->mReferenceTime));
}

template<unsigned DIM>
double VesselFlowProperties<DIM>::GetImpedanceSI(const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments) const
{
    return GetDimensionalImpedance(segments) / unit::unit_flow_impedance;
}

template<unsigned DIM>
units::quantity<unit::dynamic_viscosity> VesselFlowProperties<DIM>::GetDimensionalViscosity(const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments) const
{
    units::quantity<unit::dynamic_viscosity> value = 0.0 * unit::poiseuille;
    for (unsigned i = 0; i < segments.size(); i++)
    {
        value += segments[i]->GetFlowProperties()->GetDimensionalViscosity();
    }
    return value/ double(segments.size());
}

template<unsigned DIM>
double VesselFlowProperties<DIM>::GetViscosity(const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments) const
{
    return GetDimensionalViscosity(segments) /(this->mReferenceMass/(this->mReferenceLength * this->mReferenceTime));
}

template<unsigned DIM>
double VesselFlowProperties<DIM>::GetViscositySI(const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments) const
{
    return GetDimensionalViscosity(segments) / unit::poiseuille;
}

template<unsigned DIM>
units::quantity<unit::pressure> VesselFlowProperties<DIM>::GetDimensionalWallShearStress(const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments) const
{
    units::quantity<unit::pressure> value = 0.0 * unit::pascals;
    for (unsigned i = 0; i < segments.size(); i++)
    {
        value += segments[i]->GetFlowProperties()->GetDimensionalWallShearStress();
    }
    return value/ double(segments.size());
}

template<unsigned DIM>
double VesselFlowProperties<DIM>::GetWallShearStress(const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments) const
{
    return GetDimensionalWallShearStress(segments) /(this->mReferenceMass/(this->mReferenceLength*this->mReferenceTime*this->mReferenceTime));
}

template<unsigned DIM>
double VesselFlowProperties<DIM>::GetWallShearStressSI(const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments) const
{
    return GetDimensionalWallShearStress(segments) / unit::pascals;
}

template<unsigned DIM>
units::quantity<unit::rate> VesselFlowProperties<DIM>::GetDimensionalGrowthStimulus(const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments) const
{
    units::quantity<unit::rate> value = 0.0 * unit::reciprocal_seconds;
    for (unsigned i = 0; i < segments.size(); i++)
    {
        value += segments[i]->GetFlowProperties()->GetDimensionalGrowthStimulus();
    }
    return value/ double(segments.size());
}

template<unsigned DIM>
double VesselFlowProperties<DIM>::GetGrowthStimulus(const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments) const
{
    return GetDimensionalGrowthStimulus(segments) /(1.0/this->mReferenceTime);
}

template<unsigned DIM>
double VesselFlowProperties<DIM>::GetGrowthStimulusSI(const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments) const
{
    return GetDimensionalGrowthStimulus(segments) / unit::reciprocal_seconds;
}

template<unsigned DIM>
std::map<std::string, double> VesselFlowProperties<DIM>::GetOutputData(const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments) const
{
    std::map<std::string, double> output_data;
    output_data["Vessel Impedance kg/m^4/s"] = this->GetImpedanceSI(segments);
    output_data["Vessel Haematocrit"] = this->GetHaematocritSI(segments);
    output_data["Vessel Flow Rate m^3/s"] = this->GetFlowRateSI(segments);
    output_data["Absolute Vessel Flow Rate m^3/s"] = fabs(this->GetFlowRateSI(segments));
    output_data["Vessel Viscosity Pa.s"] = this->GetViscositySI(segments);
    output_data["Vessel Wall Shear Stress Pa"] = this->GetWallShearStressSI(segments);
    output_data["Vessel Growth Stimulus s^-1"] = this->GetGrowthStimulusSI(segments);
    return output_data;
}

template<unsigned DIM>
bool VesselFlowProperties<DIM>::HasVesselRegressed()
{
    if (SimulationTime::Instance()->GetTime()*this->mReferenceTime >= this->mRegressionTime)
    {
        assert(mUndergoingRegression);
        mRemoveViaRegression = true;
        return mRemoveViaRegression;
    }
    else
    {
        mRemoveViaRegression = false;
        return mRemoveViaRegression;
    }
}

template<unsigned DIM>
void VesselFlowProperties<DIM>::SetDimensionalHaematocrit(units::quantity<unit::dimensionless> haematocrit, const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments)
{
    for (unsigned i = 0; i < segments.size(); i++)
    {
        segments[i]->GetFlowProperties()->SetDimensionalHaematocrit(haematocrit);
    }
}

template<unsigned DIM>
void VesselFlowProperties<DIM>::SetHaematocritSI(double haematocrit, const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments)
{
    for (unsigned i = 0; i < segments.size(); i++)
    {
        segments[i]->GetFlowProperties()->SetHaematocritSI(haematocrit);
    }
}

template<unsigned DIM>
void VesselFlowProperties<DIM>::SetHaematocrit(double haematocrit, const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments)
{
    for (unsigned i = 0; i < segments.size(); i++)
    {
        segments[i]->GetFlowProperties()->SetHaematocrit(haematocrit);
    }
}

template<unsigned DIM>
void VesselFlowProperties<DIM>::SetDimensionalFlowRate(units::quantity<unit::flow_rate> haematocrit, const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments)
{
    for (unsigned i = 0; i < segments.size(); i++)
    {
        segments[i]->GetFlowProperties()->SetDimensionalFlowRate(haematocrit);
    }
}

template<unsigned DIM>
void VesselFlowProperties<DIM>::SetFlowRateSI(double haematocrit, const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments)
{
    for (unsigned i = 0; i < segments.size(); i++)
    {
        segments[i]->GetFlowProperties()->SetFlowRateSI(haematocrit);
    }
}

template<unsigned DIM>
void VesselFlowProperties<DIM>::SetFlowRate(double haematocrit, const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments)
{
    for (unsigned i = 0; i < segments.size(); i++)
    {
        segments[i]->GetFlowProperties()->SetFlowRate(haematocrit);
    }
}

template<unsigned DIM>
void VesselFlowProperties<DIM>::SetDimensionalImpedance(units::quantity<unit::flow_impedance> value, const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments)
{
    for (unsigned i = 0; i < segments.size(); i++)
    {
        segments[i]->GetFlowProperties()->SetDimensionalImpedance(value/double(segments.size()));
    }
}

template<unsigned DIM>
void VesselFlowProperties<DIM>::SetImpedanceSI(double value, const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments)
{
    for (unsigned i = 0; i < segments.size(); i++)
    {
        segments[i]->GetFlowProperties()->SetImpedanceSI(value/double(segments.size()));
    }
}

template<unsigned DIM>
void VesselFlowProperties<DIM>::SetImpedance(double value, const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments)
{
    for (unsigned i = 0; i < segments.size(); i++)
    {
        segments[i]->GetFlowProperties()->SetImpedance(value/double(segments.size()));
    }
}

template<unsigned DIM>
void VesselFlowProperties<DIM>::SetDimensionalViscosity(units::quantity<unit::dynamic_viscosity> value, const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments)
{
    for (unsigned i = 0; i < segments.size(); i++)
    {
        segments[i]->GetFlowProperties()->SetDimensionalViscosity(value);
    }
}

template<unsigned DIM>
void VesselFlowProperties<DIM>::SetViscositySI(double value, const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments)
{
    for (unsigned i = 0; i < segments.size(); i++)
    {
        segments[i]->GetFlowProperties()->SetViscositySI(value);
    }
}

template<unsigned DIM>
void VesselFlowProperties<DIM>::SetViscosity(double value, const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments)
{
    for (unsigned i = 0; i < segments.size(); i++)
    {
        segments[i]->GetFlowProperties()->SetViscosity(value);
    }
}

template<unsigned DIM>
void VesselFlowProperties<DIM>::SetDimensionalWallShearStress(units::quantity<unit::pressure> value, const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments)
{
    for (unsigned i = 0; i < segments.size(); i++)
    {
        segments[i]->GetFlowProperties()->SetDimensionalWallShearStress(value);
    }
}

template<unsigned DIM>
void VesselFlowProperties<DIM>::SetWallShearStressSI(double value, const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments)
{
    for (unsigned i = 0; i < segments.size(); i++)
    {
        segments[i]->GetFlowProperties()->SetWallShearStressSI(value);
    }
}

template<unsigned DIM>
void VesselFlowProperties<DIM>::SetWallShearStress(double value, const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments)
{
    for (unsigned i = 0; i < segments.size(); i++)
    {
        segments[i]->GetFlowProperties()->SetWallShearStress(value);
    }
}

template<unsigned DIM>
void VesselFlowProperties<DIM>::SetDimensionalGrowthStimulus(units::quantity<unit::rate> value, const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments)
{
    for (unsigned i = 0; i < segments.size(); i++)
    {
        segments[i]->GetFlowProperties()->SetDimensionalGrowthStimulus(value);
    }
}

template<unsigned DIM>
void VesselFlowProperties<DIM>::SetGrowthStimulusSI(double value, const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments)
{
    for (unsigned i = 0; i < segments.size(); i++)
    {
        segments[i]->GetFlowProperties()->SetGrowthStimulusSI(value);
    }
}

template<unsigned DIM>
void VesselFlowProperties<DIM>::SetGrowthStimulus(double value, const std::vector<boost::shared_ptr<VesselSegment<DIM> > >& segments)
{
    for (unsigned i = 0; i < segments.size(); i++)
    {
        segments[i]->GetFlowProperties()->SetGrowthStimulus(value);
    }
}

template<unsigned DIM>
void VesselFlowProperties<DIM>::SetDimensionalTimeUntilRegression(units::quantity<unit::time> time)
{
    assert(!mRemoveViaRegression);

    if (HasRegressionTimerStarted())
    {
        EXCEPTION("SetTimeUntilRegression(time) called when already undergoing regression");
    }

    mUndergoingRegression = true;
    this->mRegressionTime = SimulationTime::Instance()->GetTime()*this->mReferenceTime + time;
}

template<unsigned DIM>
void VesselFlowProperties<DIM>::SetTimeUntilRegression(double time)
{
    this->SetDimensionalTimeUntilRegression(time*this->mReferenceTime);
}

template<unsigned DIM>
bool VesselFlowProperties<DIM>::HasRegressionTimerStarted()
{
    return mUndergoingRegression;
}

template<unsigned DIM>
void VesselFlowProperties<DIM>::ResetRegressionTimer()
{
    this->mRegressionTime = DBL_MAX*this->mReferenceTime;
    mUndergoingRegression = false;
    mRemoveViaRegression = false;
}

// Explicit instantiation
template class VesselFlowProperties<2> ;
template class VesselFlowProperties<3> ;
