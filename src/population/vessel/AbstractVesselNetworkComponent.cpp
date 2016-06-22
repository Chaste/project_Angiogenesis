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

#include "AbstractVesselNetworkComponent.hpp"
#include "Exception.hpp"

template<unsigned DIM>
AbstractVesselNetworkComponent<DIM>::AbstractVesselNetworkComponent() :
        mOutputData(),
        mId(0),
        mRadius(10.0*unit::microns),
        mReferenceLength(1.0*unit::microns),
        mReferenceTime(60.0*unit::seconds),
        mReferenceMass((36.0/1e4)*unit::kg) // This leads to a 'reference pressure' of 1 Pa for the above L and T choices.
{
}

template<unsigned DIM>
AbstractVesselNetworkComponent<DIM>::~AbstractVesselNetworkComponent()
{
}

template<unsigned DIM>
unsigned AbstractVesselNetworkComponent<DIM>::GetId() const
{
    return mId;
}

template<unsigned DIM>
double AbstractVesselNetworkComponent<DIM>::GetOutputDataValue(const std::string& rKey)
{
    // Make sure the output data is up to date
    GetOutputData();

    std::map<std::string,double>::const_iterator it = mOutputData.find(rKey);
    if (it != mOutputData.end())
    {
        return it->second;
    }
    else
    {
        EXCEPTION("Requested output data key not found");
    }
}

template<unsigned DIM>
std::map<std::string, double> AbstractVesselNetworkComponent<DIM>::GetOutputData()
{
    mOutputData["VN Component Id"] = double(GetId());
    mOutputData["Dimensionless VN Component Radius"] = GetRadius();
    mOutputData["VN Component Radius: m"] = GetRadiusSI();
    return mOutputData;
}

template<unsigned DIM>
std::vector<std::string> AbstractVesselNetworkComponent<DIM>::GetOutputDataKeys()
{
    // Make sure the output data is up to date
    GetOutputData();

    std::vector<std::string> keys;
    for(std::map<std::string,double>::const_iterator it = mOutputData.begin(); it != mOutputData.end(); ++it)
    {
        keys.push_back(it->first);
    }
    return keys;
}

template<unsigned DIM>
double AbstractVesselNetworkComponent<DIM>::GetRadius() const
{
    return mRadius / mReferenceLength;
}

template<unsigned DIM>
double AbstractVesselNetworkComponent<DIM>::GetRadiusSI() const
{
    return mRadius/unit::metres;
}

template<unsigned DIM>
units::quantity<unit::length> AbstractVesselNetworkComponent<DIM>::GetDimensionalRadius() const
{
    return mRadius;
}

template<unsigned DIM>
units::quantity<unit::length> AbstractVesselNetworkComponent<DIM>::GetReferenceLength() const
{
    return mReferenceLength;
}

template<unsigned DIM>
double AbstractVesselNetworkComponent<DIM>::GetReferenceLengthSI() const
{
    return mReferenceLength/unit::metres;
}

template<unsigned DIM>
units::quantity<unit::time> AbstractVesselNetworkComponent<DIM>::GetReferenceTime() const
{
    return mReferenceTime;
}

template<unsigned DIM>
double AbstractVesselNetworkComponent<DIM>::GetReferenceTimeSI() const
{
    return mReferenceTime/unit::seconds;
}

template<unsigned DIM>
units::quantity<unit::mass> AbstractVesselNetworkComponent<DIM>::GetReferenceMass() const
{
    return mReferenceMass;
}

template<unsigned DIM>
double AbstractVesselNetworkComponent<DIM>::GetReferenceMassSI() const
{
    return mReferenceMass/unit::kg;
}

template<unsigned DIM>
void AbstractVesselNetworkComponent<DIM>::SetOutputData(const std::string& rKey, double value)
{
    mOutputData[rKey] = value;
}

template<unsigned DIM>
void AbstractVesselNetworkComponent<DIM>::SetId(unsigned id)
{
    mId = id;
}

template<unsigned DIM>
void AbstractVesselNetworkComponent<DIM>::SetDimensionalRadius(units::quantity<unit::length> radius)
{
    mRadius = radius;
}

template<unsigned DIM>
void AbstractVesselNetworkComponent<DIM>::SetRadius(double radius)
{
    mRadius = radius * mReferenceLength;
}

template<unsigned DIM>
void AbstractVesselNetworkComponent<DIM>::SetRadiusSI(double radius)
{
    mRadius = radius * unit::metres;
}

template<unsigned DIM>
void AbstractVesselNetworkComponent<DIM>::SetReferenceLength(units::quantity<unit::length> referenceLength)
{
    mReferenceLength = referenceLength;
}

template<unsigned DIM>
void AbstractVesselNetworkComponent<DIM>::SetReferenceLengthSI(double referenceLength)
{
    mReferenceLength = referenceLength * unit::metres;
}

template<unsigned DIM>
void AbstractVesselNetworkComponent<DIM>::SetReferenceTime(units::quantity<unit::time> referenceTime)
{
    mReferenceTime = referenceTime;
}

template<unsigned DIM>
void AbstractVesselNetworkComponent<DIM>::SetReferenceTimeSI(double referenceTime)
{
    mReferenceTime = referenceTime * unit::seconds;
}

template<unsigned DIM>
void AbstractVesselNetworkComponent<DIM>::SetReferenceMass(units::quantity<unit::mass> referenceMass)
{
    mReferenceMass= referenceMass;
}

template<unsigned DIM>
void AbstractVesselNetworkComponent<DIM>::SetReferenceMassSI(double referenceMass)
{
    mReferenceMass = referenceMass * unit::kg;
}

// Explicit instantiation
template class AbstractVesselNetworkComponent<2> ;
template class AbstractVesselNetworkComponent<3> ;
