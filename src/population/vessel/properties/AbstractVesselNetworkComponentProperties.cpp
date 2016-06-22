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

#include "AbstractVesselNetworkComponentProperties.hpp"

template<unsigned DIM>
AbstractVesselNetworkComponentProperties<DIM>::AbstractVesselNetworkComponentProperties() :
        mReferenceLength(1.0*unit::microns),
        mReferenceTime(60.0*unit::seconds),
        mReferenceMass((36.0/1e4)*unit::kg) // This leads to a 'reference pressure' of 1 Pa for the above L and T choices, which is nice.
{
}

template<unsigned DIM>
AbstractVesselNetworkComponentProperties<DIM>::~AbstractVesselNetworkComponentProperties()
{
}

template<unsigned DIM>
std::map<std::string, double> AbstractVesselNetworkComponentProperties<DIM>::GetOutputData() const
{
    std::map<std::string, double> output_data;
    output_data["Reference Length m"] = mReferenceLength/unit::metres;
    output_data["Reference Mass s"] = mReferenceTime/unit::seconds;
    output_data["Reference Time kg"] = mReferenceMass/unit::kg;
    return output_data;
}

template<unsigned DIM>
units::quantity<unit::length> AbstractVesselNetworkComponentProperties<DIM>::GetReferenceLength() const
{
    return mReferenceLength;
}

template<unsigned DIM>
units::quantity<unit::time> AbstractVesselNetworkComponentProperties<DIM>::GetReferenceTime() const
{
    return mReferenceTime;
}

template<unsigned DIM>
units::quantity<unit::mass> AbstractVesselNetworkComponentProperties<DIM>::GetReferenceMass() const
{
    return mReferenceMass;
}

template<unsigned DIM>
double AbstractVesselNetworkComponentProperties<DIM>::GetReferenceLengthSI() const
{
    return mReferenceLength/unit::metres;
}

template<unsigned DIM>
double AbstractVesselNetworkComponentProperties<DIM>::GetReferenceMassSI() const
{
    return mReferenceMass/unit::kg;
}

template<unsigned DIM>
double AbstractVesselNetworkComponentProperties<DIM>::GetReferenceTimeSI() const
{
    return mReferenceTime/unit::seconds;
}

template<unsigned DIM>
void AbstractVesselNetworkComponentProperties<DIM>::SetReferenceLength(units::quantity<unit::length> referenceLength)
{
    mReferenceLength = referenceLength;
}

template<unsigned DIM>
void AbstractVesselNetworkComponentProperties<DIM>::SetReferenceTime(units::quantity<unit::time> referenceTime)
{
    mReferenceTime = referenceTime;
}

template<unsigned DIM>
void AbstractVesselNetworkComponentProperties<DIM>::SetReferenceMass(units::quantity<unit::mass> referenceMass)
{
    mReferenceMass = referenceMass;
}

template<unsigned DIM>
void AbstractVesselNetworkComponentProperties<DIM>::SetReferenceLengthSI(double referenceLength)
{
    mReferenceLength = referenceLength * unit::metres;
}

template<unsigned DIM>
void AbstractVesselNetworkComponentProperties<DIM>::SetReferenceMassSI(double referenceMass)
{
    mReferenceMass = referenceMass * unit::kg;
}

template<unsigned DIM>
void AbstractVesselNetworkComponentProperties<DIM>::SetReferenceTimeSI(double referenceTime)
{
    mReferenceTime= referenceTime * unit::seconds;
}

// Explicit instantiation
template class AbstractVesselNetworkComponentProperties<2>;
template class AbstractVesselNetworkComponentProperties<3>;
