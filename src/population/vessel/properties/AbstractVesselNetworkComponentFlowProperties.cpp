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

#include "AbstractVesselNetworkComponentFlowProperties.hpp"

template<unsigned DIM>
AbstractVesselNetworkComponentFlowProperties<DIM>::AbstractVesselNetworkComponentFlowProperties() : AbstractVesselNetworkComponentProperties<DIM>(),
        mPressure(0.0 * unit::pascals)
{
}

template<unsigned DIM>
AbstractVesselNetworkComponentFlowProperties<DIM>::~AbstractVesselNetworkComponentFlowProperties()
{
}

template<unsigned DIM>
double AbstractVesselNetworkComponentFlowProperties<DIM>::GetPressure() const
{
    return mPressure/(this->mReferenceMass/(this->mReferenceLength*this->mReferenceTime*this->mReferenceTime));
}

template<unsigned DIM>
units::quantity<unit::pressure> AbstractVesselNetworkComponentFlowProperties<DIM>::GetDimensionalPressure() const
{
    return mPressure;
}

template<unsigned DIM>
double AbstractVesselNetworkComponentFlowProperties<DIM>::GetPressureSI() const
{
    return mPressure/unit::pascals;
}

template<unsigned DIM>
std::map<std::string, double> AbstractVesselNetworkComponentFlowProperties<DIM>::GetOutputData() const
{
    std::map<std::string, double> output_data;
    output_data["Reference Length m"] = this->GetReferenceLengthSI();
    output_data["Reference Mass kg"] = this->GetReferenceMassSI();
    output_data["Reference Time s"] = this->GetReferenceTimeSI();
    return output_data;
}

template<unsigned DIM>
void AbstractVesselNetworkComponentFlowProperties<DIM>::SetPressure(double pressure)
{
    mPressure = pressure * (this->mReferenceMass/(this->mReferenceLength*this->mReferenceTime*this->mReferenceTime));
}

template<unsigned DIM>
void AbstractVesselNetworkComponentFlowProperties<DIM>::SetDimensionalPressure(units::quantity<unit::pressure> pressure)
{
    mPressure = pressure;
}

template<unsigned DIM>
void AbstractVesselNetworkComponentFlowProperties<DIM>::SetPressureSI(double pressure)
{
    mPressure = pressure * unit::pascals;
}

// Explicit instantiation
template class AbstractVesselNetworkComponentFlowProperties<2> ;
template class AbstractVesselNetworkComponentFlowProperties<3> ;
