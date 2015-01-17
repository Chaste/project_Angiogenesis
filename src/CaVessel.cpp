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

#include "CaVessel.hpp"

template<unsigned DIM>
CaVessel<DIM>::CaVessel()
	: pNode1(new CaVascularNetworkNode<DIM>()),
	  pNode2(new CaVascularNetworkNode<DIM>()),
	  mActiveTipCellAtNode1(false),
	  mActiveTipCellAtNode2(false),
	  mVesselSegmentLocations(),
	  mTimeWithLowWallShearStress(0.0),
	  mRadius(0.0),
	  mPreviousRadius(0.0),
	  mHaematocritLevel(0.0),
	  mFlowVelocity(0.0),
	  mFlowRate(0.0),
	  mImpedance(0.0),
	  mLength(0.0),
	  mWallShearStress(0.0),
	  mViscosity(0.0),
	  mMechanicalStimulus(0.0),
	  mMetabolicStimulus(0.0),
	  mShrinkingStimulus(0.0),
	  mDownstreamConductedStimulus(0.0),
	  mUpstreamConductedStimulus(0.0),
	  mChemicalCollection(),
	  mIsPartOfNeovasculature(true),
	  mCanExtend(false)
{
}

template<unsigned DIM>
CaVessel<DIM>::~CaVessel()
{
}

template<unsigned DIM>
bool CaVessel<DIM>::IsInputVessel()
{
    return (pNode1->IsInputNode() || pNode2->IsInputNode());
}


template<unsigned DIM>
void CaVessel<DIM>::CopyMechanicalPropertyValuesAndChemicalConcentrations(boost::shared_ptr<CaVessel<DIM> > another_vessel)
{
    mTimeWithLowWallShearStress = another_vessel->GetTimeWithLowWallShearStress();
    mRadius = another_vessel->GetRadius();
    mPreviousRadius = another_vessel->GetPreviousRadius();
    mHaematocritLevel = another_vessel->GetHaematocritLevel();
    mFlowVelocity = another_vessel->GetFlowVelocity();
    mFlowRate = another_vessel->GetFlowRate();
    mImpedance = another_vessel->GetImpedance();
    mLength = another_vessel->GetLength();
    mWallShearStress = another_vessel->GetWallShearStress();
    mViscosity = another_vessel->GetViscosity();
    mMechanicalStimulus = another_vessel->GetMechanicalStimulus();
    mMetabolicStimulus = another_vessel->GetMetabolicStimulus();
    mShrinkingStimulus = another_vessel->GetShrinkingStimulus();
    mDownstreamConductedStimulus = another_vessel->GetDownstreamConductedStimulus();
    mUpstreamConductedStimulus = another_vessel->GetUpstreamConductedStimulus();
    mCanExtend = another_vessel->CanExtend();

    for (unsigned i = 0; i < another_vessel->GetNumberOfIntraVascularChemicals(); i++)
    {
        mChemicalCollection.AddIntraVascularChemical(another_vessel->rGetCollectionOfIntraVascularChemicals().GetIntraVascularChemicalCollection()[i].GetChemicalName(),
                Concentration(another_vessel->rGetCollectionOfIntraVascularChemicals().GetIntraVascularChemicalCollection()[i].GetConcentration(),
                another_vessel->rGetCollectionOfIntraVascularChemicals().GetIntraVascularChemicalCollection()[i].GetUnits()),
                another_vessel->rGetCollectionOfIntraVascularChemicals().GetIntraVascularChemicalCollection()[i].GetPermeabilityOfVesselWallToChemical());
    }
}

template<unsigned DIM>
double CaVessel<DIM>::GetRadius()
{
    return mRadius;
}

template<unsigned DIM>
double CaVessel<DIM>::GetTortuosity()
{
    double distance_between_ends_of_vessel;

    if (DIM == 3u)
    {
        distance_between_ends_of_vessel = pow((pow(double(GetNode1()->GetLocation()[0] - GetNode2()->GetLocation()[0]), 2) +
                pow(double(GetNode1()->GetLocation()[1] - GetNode2()->GetLocation()[1]), 2) +
                pow(double(GetNode1()->GetLocation()[2] - GetNode2()->GetLocation()[2]), 2)), 0.5);
    }
    else
    {
		distance_between_ends_of_vessel = pow((pow(double(GetNode1()->GetLocation()[0] - GetNode2()->GetLocation()[0]), 2) +
				pow(double(GetNode1()->GetLocation()[1] - GetNode2()->GetLocation()[1]), 2)), 0.5);
    }

    return (mLength/distance_between_ends_of_vessel);
}

template<unsigned DIM>
double CaVessel<DIM>::GetPreviousRadius()
{
    return mPreviousRadius;
}

template<unsigned DIM>
double CaVessel<DIM>::GetHaematocritLevel()
{
    return mHaematocritLevel;
}

template<unsigned DIM>
double CaVessel<DIM>::GetFlowVelocity()
{
    return mFlowVelocity;
}

template<unsigned DIM>
double CaVessel<DIM>::GetFlowRate()
{
    return mFlowRate;
}

template<unsigned DIM>
double CaVessel<DIM>::GetImpedance()
{
    return mImpedance;
}

template<unsigned DIM>
double CaVessel<DIM>::GetLength()
{
    return mLength;
}

template<unsigned DIM>
double CaVessel<DIM>::GetWallShearStress()
{
    return mWallShearStress;
}

template<unsigned DIM>
double CaVessel<DIM>::GetViscosity()
{
    return mViscosity;
}

template<unsigned DIM>
double CaVessel<DIM>::GetMechanicalStimulus()
{
    return mMechanicalStimulus;
}

template<unsigned DIM>
double CaVessel<DIM>::GetMetabolicStimulus()
{
    return mMetabolicStimulus;
}

template<unsigned DIM>
double CaVessel<DIM>::GetShrinkingStimulus()
{
    return mShrinkingStimulus;
}

template<unsigned DIM>
double CaVessel<DIM>::GetDownstreamConductedStimulus()
{
    return mDownstreamConductedStimulus;
}

template<unsigned DIM>
double CaVessel<DIM>::GetUpstreamConductedStimulus()
{
    return mUpstreamConductedStimulus;
}

template<unsigned DIM>
boost::shared_ptr<CaVessel<DIM> > CaVessel<DIM>::shared()
{
    return this->shared_from_this();
}

template<unsigned DIM>
boost::shared_ptr<CaVascularNetworkNode<DIM> > CaVessel<DIM>::GetNode1()
{
    return pNode1;
}

template<unsigned DIM>
boost::shared_ptr<CaVascularNetworkNode<DIM> > CaVessel<DIM>::GetNode2()
{
    return pNode2;
}

template<unsigned DIM>
double CaVessel<DIM>::GetTimeWithLowWallShearStress()
{
    return mTimeWithLowWallShearStress;
}

template<unsigned DIM>
std::vector<ChastePoint<DIM> > CaVessel<DIM>::GetSegmentCoordinates()
{
    return mVesselSegmentLocations;
}

template<unsigned DIM>
ChastePoint<DIM> CaVessel<DIM>::GetSegmentCoordinate(unsigned i)
{
	if(i < mVesselSegmentLocations.size())
	{
		return mVesselSegmentLocations[i];
	}
	else
	{
		EXCEPTION("Attempted to access a vessel segment with an out of range index.");
	}
}

template<unsigned DIM>
unsigned CaVessel<DIM>::GetNumberOfSegments()
{
    return mVesselSegmentLocations.size();
}


template<unsigned DIM>
bool CaVessel<DIM>::HasActiveTipCell()
{
    return (ActiveTipCellLocatedAtNode1() || ActiveTipCellLocatedAtNode2());
}

template<unsigned DIM>
bool CaVessel<DIM>::ActiveTipCellLocatedAtNode1()
{
    return mActiveTipCellAtNode1;
}

template<unsigned DIM>
bool CaVessel<DIM>::ActiveTipCellLocatedAtNode2()
{
    return mActiveTipCellAtNode2;
}

template<unsigned DIM>
IntraVascularChemicalCollection& CaVessel<DIM>::rGetCollectionOfIntraVascularChemicals()
{
    return mChemicalCollection;
}

template<unsigned DIM>
unsigned CaVessel<DIM>::GetNumberOfIntraVascularChemicals()
{
    return mChemicalCollection.GetIntraVascularChemicalCollection().size();
}

template<unsigned DIM>
double CaVessel<DIM>::GetIntraVascularChemicalConcentration(string chemical_name)
{
    return mChemicalCollection.GetIntraVascularChemicalConcentration(chemical_name);
}


template<unsigned DIM>
bool CaVessel<DIM>::IsPartOfNeovasculature()
{
    return mIsPartOfNeovasculature;
}

template<unsigned DIM>
bool CaVessel<DIM>::CanExtend()
{
    return mCanExtend;
}


template<unsigned DIM>
void CaVessel<DIM>::SetIntraVascularChemicalConcentration(string chemical_name, Concentration concentration)
{
    mChemicalCollection.SetIntraVascularChemicalConcentration(chemical_name, concentration);
}

template<unsigned DIM>
void CaVessel<DIM>::SetNode1(boost::shared_ptr<CaVascularNetworkNode<DIM> > node)
{
    pNode1 = node;
}

template<unsigned DIM>
void CaVessel<DIM>::SetNode1Location(ChastePoint<DIM> location)
{
    pNode1->SetLocation(location);
    pNode1->AddAdjoiningVessel(shared()); // badly placed - doesn't make much sense but only place to do it really
}

template<unsigned DIM>
void CaVessel<DIM>::SetNode2(boost::shared_ptr<CaVascularNetworkNode<DIM> > node)
{
    pNode2 = node;
}

template<unsigned DIM>
void CaVessel<DIM>::SetNode2Location(ChastePoint<DIM> location)
{
    pNode2->SetLocation(location);
    pNode2->AddAdjoiningVessel(shared());
}

/*
 * todo This method needs some adjustment ... dealing with coordinates that may not be integers, for example.
 * And what if there is no z-coordinate
 */
template<unsigned DIM>
void CaVessel<DIM>::SetNextVesselSegmentCoordinate(ChastePoint<DIM> point)
{
    if (GetNumberOfSegments() > 0)
    {
        mVesselSegmentLocations.push_back(point);
        SetLength(0.0);

        for (unsigned i = 1; i < GetNumberOfSegments(); i++)
        {
            double segment_length = 0;

            if (DIM == 2)
            {
                segment_length = pow((pow(double(GetSegmentCoordinate(i)[0]-GetSegmentCoordinate(i-1)[0]), 2) +
                		pow(double(GetSegmentCoordinate(i)[1]-GetSegmentCoordinate(i-1)[1]), 2)), 0.5);
            }

            if (DIM == 3)
            {
                segment_length = pow((pow(double(GetSegmentCoordinate(i)[0]-GetSegmentCoordinate(i-1)[0]), 2) +
                		pow(double(GetSegmentCoordinate(i)[1]-GetSegmentCoordinate(i-1)[1]),2) +
                		pow(double(GetSegmentCoordinate(i)[2] - GetSegmentCoordinate(i - 1)[2]), 2)), 0.5);
            }
            SetLength(GetLength()+segment_length);
        }
    }
    else
    {
        mVesselSegmentLocations.push_back(point);
        SetLength(1.0);
    }
}

template<unsigned DIM>
void CaVessel<DIM>::SetTimeWithLowWallShearStress(double time)
{
    mTimeWithLowWallShearStress = time;
}

template<unsigned DIM>
void CaVessel<DIM>::IncrementTimeWithLowWallShearStress(double delta_time)
{
    mTimeWithLowWallShearStress += delta_time;
}

template<unsigned DIM>
void CaVessel<DIM>::SetActiveTipCellLocatedAtNode1(bool value)
{
    if (value == true)
    {
    	if (!(GetNode1()->GetNumberOfAdjoiningVessels() == 1 && GetNode1()->GetAdjoiningVessel(0) == shared()))
    	{
    		EXCEPTION("Node 1 is not linked with the vessel or there is an incorrect number of vessels attached to the node.");
    	}

    	if (mActiveTipCellAtNode2)
    	{
    		EXCEPTION("There is already an active tip cell at the other end of the vessel.");
    	}
    }
    mActiveTipCellAtNode1 = value;
}

template<unsigned DIM>
void CaVessel<DIM>::SetActiveTipCellLocatedAtNode2(bool value)
{
    if (value == true)
    {

    	if (!(GetNode2()->GetNumberOfAdjoiningVessels() == 1 && GetNode2()->GetAdjoiningVessel(0) == shared()))
    	{
    		EXCEPTION("Node 2 is not linked with the vessel or there is an incorrect number of vessels attached to the node.");
    	}

    	if (mActiveTipCellAtNode1)
    	{
    		EXCEPTION("There is already an active tip cell at the other end of the vessel.");
    	}
    }
    mActiveTipCellAtNode2 = value;
}

template<unsigned DIM>
void CaVessel<DIM>::SetRadius(double value)
{
    mRadius = value;
}

template<unsigned DIM>
void CaVessel<DIM>::SetPreviousRadius(double value)
{
    mPreviousRadius = value;
}

template<unsigned DIM>
void CaVessel<DIM>::SetHaematocritLevel(double value)
{
    mHaematocritLevel = value;
}

template<unsigned DIM>
void CaVessel<DIM>::SetFlowVelocity(double value)
{
    mFlowVelocity = value;
}

template<unsigned DIM>
void CaVessel<DIM>::SetFlowRate(double value)
{
    mFlowRate = value;
}

template<unsigned DIM>
void CaVessel<DIM>::SetImpedance(double value)
{
    mImpedance = value;
}

template<unsigned DIM>
void CaVessel<DIM>::SetLength(double value)
{
    mLength = value;
}

template<unsigned DIM>
void CaVessel<DIM>::SetWallShearStress(double value)
{
    mWallShearStress = value;
}

template<unsigned DIM>
void CaVessel<DIM>::SetViscosity(double value)
{
    mViscosity = value;
}

template<unsigned DIM>
void CaVessel<DIM>::SetMechanicalStimulus(double value)
{
    mMechanicalStimulus = value;
}

template<unsigned DIM>
void CaVessel<DIM>::SetMetabolicStimulus(double value)
{
    mMetabolicStimulus = value;
}

template<unsigned DIM>
void CaVessel<DIM>::SetShrinkingStimulus(double value)
{
    mShrinkingStimulus = value;
}

template<unsigned DIM>
void CaVessel<DIM>::SetDownstreamConductedStimulus(double value)
{
    mDownstreamConductedStimulus = value;
}

template<unsigned DIM>
void CaVessel<DIM>::SetUpstreamConductedStimulus(double value)
{
    mUpstreamConductedStimulus = value;
}

template<unsigned DIM>
bool CaVessel<DIM>::IsAttachedToNode(boost::shared_ptr<CaVascularNetworkNode<DIM> > node)
{
    return (node == pNode1 || node == pNode2);
}

template<unsigned DIM>
void CaVessel<DIM>::SetIsPartOfNeovasculature(bool value)
{
    mIsPartOfNeovasculature = value;
}

template<unsigned DIM>
void CaVessel<DIM>::CanExtend(bool value)
{
    mCanExtend = value;
}

// Explicit instantiation
template class CaVessel<2>;
template class CaVessel<3>;
