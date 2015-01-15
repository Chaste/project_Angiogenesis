/*
 * CaVessel.cpp
 *
 *  Created on: 13 Jan 2015
 *      Author: connor
 */

#include "CaVessel.hpp"

template<unsigned SPATIAL_DIM>
CaVessel<SPATIAL_DIM>::CaVessel() : pNode1(new CaVascularNetworkNode<SPATIAL_DIM>()),pNode2(new CaVascularNetworkNode<SPATIAL_DIM>()),mActiveTipCellLocatedAtNode1(false),mActiveTipCellLocatedAtNode2(false),mVesselSegmentLocations(),mTimeWithLowWallShearStress(0.0),mRadius(1.0*pow(10.0,-5)),mPreviousRadius(1.0),mHaematocritLevel(0.45),mFlowVelocity(0.0),mFlowRate(0.0),mImpedance(0.0),mLength(0.0),mWallShearStress(0.0),mViscosity(0.0),mMechanicalStimulus(0.0),mMetabolicStimulus(0.0),mShrinkingStimulus(0.0),mDownstreamConductedStimulus(0.0),mUpstreamConductedStimulus(0.0),mChemicalCollection(),mIsPartOfNeovasculature(true), mCanExtend(false)
{

}

template<unsigned SPATIAL_DIM>
CaVessel<SPATIAL_DIM>::~CaVessel()
{

}

template<unsigned SPATIAL_DIM>
bool CaVessel<SPATIAL_DIM>::IsInputVessel()
{
    return (pNode1->IsInputNode() || pNode2->IsInputNode());
}


template<unsigned SPATIAL_DIM>
void CaVessel<SPATIAL_DIM>::CopyMechanicalPropertyValuesAndChemicalConcentrations(boost::shared_ptr<CaVessel<SPATIAL_DIM> > anotherVessel)
{
    mTimeWithLowWallShearStress = anotherVessel->GetTimeWithLowWallShearStress();
    mRadius = anotherVessel->GetRadius();
    mPreviousRadius = anotherVessel->GetPreviousRadius();
    mHaematocritLevel = anotherVessel->GetHaematocritLevel();
    mFlowVelocity = anotherVessel->GetFlowVelocity();
    mFlowRate = anotherVessel->GetFlowRate();
    mImpedance = anotherVessel->GetImpedance();
    mLength = anotherVessel->GetLength();
    mWallShearStress = anotherVessel->GetWallShearStress();
    mViscosity = anotherVessel->GetViscosity();
    mMechanicalStimulus = anotherVessel->GetMechanicalStimulus();
    mMetabolicStimulus = anotherVessel->GetMetabolicStimulus();
    mShrinkingStimulus = anotherVessel->GetShrinkingStimulus();
    mDownstreamConductedStimulus = anotherVessel->GetDownstreamConductedStimulus();
    mUpstreamConductedStimulus = anotherVessel->GetUpstreamConductedStimulus();
    mCanExtend = anotherVessel->CanExtend();

    for (unsigned i = 0; i < anotherVessel->GetNumberOfIntraVascularChemicals(); i++)
    {
        mChemicalCollection.AddIntraVascularChemical(anotherVessel->GetCollectionOfIntraVascularChemicals().GetIntraVascularChemicalCollection()[i].GetChemicalName(),
                                                     Concentration(anotherVessel->GetCollectionOfIntraVascularChemicals().GetIntraVascularChemicalCollection()[i].GetConcentration(),
                                                                   anotherVessel->GetCollectionOfIntraVascularChemicals().GetIntraVascularChemicalCollection()[i].GetUnits()),
                                                                   anotherVessel->GetCollectionOfIntraVascularChemicals().GetIntraVascularChemicalCollection()[i].GetPermeabilityOfVesselWallToChemical());
    }

}

template<unsigned SPATIAL_DIM>
double CaVessel<SPATIAL_DIM>::GetRadius()
{
    return mRadius;
}

template<unsigned SPATIAL_DIM>
double CaVessel<SPATIAL_DIM>::GetTortuosity()
{
    double distanceBetweenEndsOfVessel;

    if (SPATIAL_DIM > 2)
    {
        distanceBetweenEndsOfVessel = pow((pow(double(GetNode1()->GetLocation()[0] - GetNode2()->GetLocation()[0]),2) +
                pow(double(GetNode1()->GetLocation()[1] - GetNode2()->GetLocation()[1]),2) +
                pow(double(GetNode1()->GetLocation()[2] - GetNode2()->GetLocation()[2]),2)),0.5);
    }
    else
        {
            distanceBetweenEndsOfVessel = pow((pow(double(GetNode1()->GetLocation()[0] - GetNode2()->GetLocation()[0]),2) +
                    pow(double(GetNode1()->GetLocation()[1] - GetNode2()->GetLocation()[1]),2)),0.5);
        }

    return (mLength/distanceBetweenEndsOfVessel);
}

template<unsigned SPATIAL_DIM>
double CaVessel<SPATIAL_DIM>::GetPreviousRadius()
{
    return mPreviousRadius;
}

template<unsigned SPATIAL_DIM>
double CaVessel<SPATIAL_DIM>::GetHaematocritLevel()
{
    return mHaematocritLevel;
}

template<unsigned SPATIAL_DIM>
double CaVessel<SPATIAL_DIM>::GetFlowVelocity()
{
    return mFlowVelocity;
}

template<unsigned SPATIAL_DIM>
double CaVessel<SPATIAL_DIM>::GetFlowRate()
{
    return mFlowRate;
}

template<unsigned SPATIAL_DIM>
double CaVessel<SPATIAL_DIM>::GetImpedance()
{
    return mImpedance;
}

template<unsigned SPATIAL_DIM>
double CaVessel<SPATIAL_DIM>::GetLength()
{
    return mLength;
}

template<unsigned SPATIAL_DIM>
double CaVessel<SPATIAL_DIM>::GetWallShearStress()
{
    return mWallShearStress;
}

template<unsigned SPATIAL_DIM>
double CaVessel<SPATIAL_DIM>::GetViscosity()
{
    return mViscosity;
}

template<unsigned SPATIAL_DIM>
double CaVessel<SPATIAL_DIM>::GetMechanicalStimulus()
{
    return mMechanicalStimulus;
}

template<unsigned SPATIAL_DIM>
double CaVessel<SPATIAL_DIM>::GetMetabolicStimulus()
{
    return mMetabolicStimulus;
}

template<unsigned SPATIAL_DIM>
double CaVessel<SPATIAL_DIM>::GetShrinkingStimulus()
{
    return mShrinkingStimulus;
}

template<unsigned SPATIAL_DIM>
double CaVessel<SPATIAL_DIM>::GetDownstreamConductedStimulus()
{
    return mDownstreamConductedStimulus;
}

template<unsigned SPATIAL_DIM>
double CaVessel<SPATIAL_DIM>::GetUpstreamConductedStimulus()
{
    return mUpstreamConductedStimulus;
}

template<unsigned SPATIAL_DIM>
boost::shared_ptr<CaVessel<SPATIAL_DIM> > CaVessel<SPATIAL_DIM>::shared()
{
    return this->shared_from_this();
}

template<unsigned SPATIAL_DIM>
boost::shared_ptr<CaVascularNetworkNode<SPATIAL_DIM> > CaVessel<SPATIAL_DIM>::GetNode1()
{
    return pNode1;
}

template<unsigned SPATIAL_DIM>
boost::shared_ptr<CaVascularNetworkNode<SPATIAL_DIM> > CaVessel<SPATIAL_DIM>::GetNode2()
{
    return pNode2;
}

template<unsigned SPATIAL_DIM>
double CaVessel<SPATIAL_DIM>::GetTimeWithLowWallShearStress() const
{
    return mTimeWithLowWallShearStress;
}

template<unsigned SPATIAL_DIM>
std::vector<ChastePoint<SPATIAL_DIM> > CaVessel<SPATIAL_DIM>::GetSegmentCoordinates()
{
    return mVesselSegmentLocations;
}

template<unsigned SPATIAL_DIM>
ChastePoint<SPATIAL_DIM> CaVessel<SPATIAL_DIM>::GetSegmentCoordinate(unsigned i)
{
    assert(i < mVesselSegmentLocations.size());
    return mVesselSegmentLocations[i];
}

template<unsigned SPATIAL_DIM>
unsigned CaVessel<SPATIAL_DIM>::GetNumberOfSegments()
{
    return mVesselSegmentLocations.size();
}


template<unsigned SPATIAL_DIM>
bool CaVessel<SPATIAL_DIM>::HasActiveTipCell() const
{
    return (ActiveTipCellLocatedAtNode1() || ActiveTipCellLocatedAtNode2());
}

template<unsigned SPATIAL_DIM>
bool CaVessel<SPATIAL_DIM>::ActiveTipCellLocatedAtNode1() const
{
    return mActiveTipCellLocatedAtNode1;
}

template<unsigned SPATIAL_DIM>
bool CaVessel<SPATIAL_DIM>::ActiveTipCellLocatedAtNode2() const
{
    return mActiveTipCellLocatedAtNode2;
}

template<unsigned SPATIAL_DIM>
IntraVascularChemicalCollection& CaVessel<SPATIAL_DIM>::GetCollectionOfIntraVascularChemicals()
{
    return mChemicalCollection;
}

template<unsigned SPATIAL_DIM>
unsigned CaVessel<SPATIAL_DIM>::GetNumberOfIntraVascularChemicals()
{
    return mChemicalCollection.GetIntraVascularChemicalCollection().size();
}

template<unsigned SPATIAL_DIM>
double CaVessel<SPATIAL_DIM>::GetIntraVascularChemicalConcentration(string chemicalname)
{
    return mChemicalCollection.GetIntraVascularChemicalConcentration(chemicalname);
}


template<unsigned SPATIAL_DIM>
bool CaVessel<SPATIAL_DIM>::IsPartOfNeovasculature()
{
    return mIsPartOfNeovasculature;
}

template<unsigned SPATIAL_DIM>
bool CaVessel<SPATIAL_DIM>::CanExtend()
{
    return mCanExtend;
}


template<unsigned SPATIAL_DIM>
void CaVessel<SPATIAL_DIM>::SetIntraVascularChemicalConcentration(string chemicalname, Concentration concentration)
{
    mChemicalCollection.SetIntraVascularChemicalConcentration(chemicalname,concentration);
}

template<unsigned SPATIAL_DIM>
void CaVessel<SPATIAL_DIM>::SetNode1(boost::shared_ptr<CaVascularNetworkNode<SPATIAL_DIM> > node)
{
    pNode1 = node;
}

template<unsigned SPATIAL_DIM>
void CaVessel<SPATIAL_DIM>::SetNode1Location(ChastePoint<SPATIAL_DIM> location)
{
    pNode1->SetLocation(location);
    pNode1->AddAdjoiningVessel(shared()); // badly placed - doesn't make much sense but only place to do it really
}

template<unsigned SPATIAL_DIM>
void CaVessel<SPATIAL_DIM>::SetNode2(boost::shared_ptr<CaVascularNetworkNode<SPATIAL_DIM> > node)
{
    pNode2 = node;
}

template<unsigned SPATIAL_DIM>
void CaVessel<SPATIAL_DIM>::SetNode2Location(ChastePoint<SPATIAL_DIM> location)
{
    pNode2->SetLocation(location);
    pNode2->AddAdjoiningVessel(shared());
}

/*
 * todo This method needs some adjustment ... dealing with coordinates that may not be integers, for example. And what if there is no z-coordinate
 */
template<unsigned SPATIAL_DIM>
void CaVessel<SPATIAL_DIM>::SetNextVesselSegmentCoordinate(ChastePoint<SPATIAL_DIM> Coords)
{

    //    // method also includes re-calculation of length of vessel, given that vessel has had another segment added
    //
    //    // commented out code would enable node locations to be automatically updated - if implementing this need to change how vessel networks are constructed and how vessel networks mutate (i.e. how vessels are extended, etc.
    //
    if (GetNumberOfSegments() > 0)
    {
        //
        //        // check that coordinate of next vessel segment is adjacent to last segment coordinates.
        //
        //        bool adjacent = false;
        //
        //        if (Coords[0] == GetSegmentCoordinate(GetNumberOfSegments() - 1)[0] + 1 && Coords[1] == GetSegmentCoordinate(GetNumberOfSegments() - 1)[1] && Coords[2] == GetSegmentCoordinate(GetNumberOfSegments() - 1)[2])
        //        {
        //            adjacent = true;
        //        }
        //        if (Coords[0] == GetSegmentCoordinate(GetNumberOfSegments() - 1)[0] - 1 && Coords[1] == GetSegmentCoordinate(GetNumberOfSegments() - 1)[1] && Coords[2] == GetSegmentCoordinate(GetNumberOfSegments() - 1)[2])
        //        {
        //            adjacent = true;
        //        }
        //        if (Coords[0] == GetSegmentCoordinate(GetNumberOfSegments() - 1)[0] + 1 && Coords[1] == GetSegmentCoordinate(GetNumberOfSegments() - 1)[1] + 1 && Coords[2] == GetSegmentCoordinate(GetNumberOfSegments() - 1)[2])
        //        {
        //            adjacent = true;
        //        }
        //        if (Coords[0] == GetSegmentCoordinate(GetNumberOfSegments() - 1)[0] + 1 && Coords[1] == GetSegmentCoordinate(GetNumberOfSegments() - 1)[1] - 1 && Coords[2] == GetSegmentCoordinate(GetNumberOfSegments() - 1)[2])
        //        {
        //            adjacent = true;
        //        }
        //        if (Coords[0] == GetSegmentCoordinate(GetNumberOfSegments() - 1)[0] - 1 && Coords[1] == GetSegmentCoordinate(GetNumberOfSegments() - 1)[1] + 1 && Coords[2] == GetSegmentCoordinate(GetNumberOfSegments() - 1)[2])
        //        {
        //            adjacent = true;
        //        }
        //        if (Coords[0] == GetSegmentCoordinate(GetNumberOfSegments() - 1)[0] - 1 && Coords[1] == GetSegmentCoordinate(GetNumberOfSegments() - 1)[1] - 1 && Coords[2] == GetSegmentCoordinate(GetNumberOfSegments() - 1)[2])
        //        {
        //            adjacent = true;
        //        }
        //        if (Coords[0] == GetSegmentCoordinate(GetNumberOfSegments() - 1)[0] + 1 && Coords[1] == GetSegmentCoordinate(GetNumberOfSegments() - 1)[1] - 1 && Coords[2] == GetSegmentCoordinate(GetNumberOfSegments() - 1)[2] - 1)
        //        {
        //            adjacent = true;
        //        }
        //        if (Coords[0] == GetSegmentCoordinate(GetNumberOfSegments() - 1)[0] + 1 && Coords[1] == GetSegmentCoordinate(GetNumberOfSegments() - 1)[1] - 1 && Coords[2] == GetSegmentCoordinate(GetNumberOfSegments() - 1)[2] + 1)
        //        {
        //            adjacent = true;
        //        }
        //        if (Coords[0] == GetSegmentCoordinate(GetNumberOfSegments() - 1)[0] + 1 && Coords[1] == GetSegmentCoordinate(GetNumberOfSegments() - 1)[1] + 1 && Coords[2] == GetSegmentCoordinate(GetNumberOfSegments() - 1)[2] - 1)
        //        {
        //            adjacent = true;
        //        }
        //        if (Coords[0] == GetSegmentCoordinate(GetNumberOfSegments() - 1)[0] + 1 && Coords[1] == GetSegmentCoordinate(GetNumberOfSegments() - 1)[1] + 1 && Coords[2] == GetSegmentCoordinate(GetNumberOfSegments() - 1)[2] + 1)
        //        {
        //            adjacent = true;
        //        }
        //        if (Coords[0] == GetSegmentCoordinate(GetNumberOfSegments() - 1)[0] - 1 && Coords[1] == GetSegmentCoordinate(GetNumberOfSegments() - 1)[1] - 1 && Coords[2] == GetSegmentCoordinate(GetNumberOfSegments() - 1)[2] - 1)
        //        {
        //            adjacent = true;
        //        }
        //        if (Coords[0] == GetSegmentCoordinate(GetNumberOfSegments() - 1)[0] - 1 && Coords[1] == GetSegmentCoordinate(GetNumberOfSegments() - 1)[1] - 1 && Coords[2] == GetSegmentCoordinate(GetNumberOfSegments() - 1)[2] + 1)
        //        {
        //            adjacent = true;
        //        }
        //        if (Coords[0] == GetSegmentCoordinate(GetNumberOfSegments() - 1)[0] - 1 && Coords[1] == GetSegmentCoordinate(GetNumberOfSegments() - 1)[1] + 1 && Coords[2] == GetSegmentCoordinate(GetNumberOfSegments() - 1)[2] - 1)
        //        {
        //            adjacent = true;
        //        }
        //        if (Coords[0] == GetSegmentCoordinate(GetNumberOfSegments() - 1)[0] - 1 && Coords[1] == GetSegmentCoordinate(GetNumberOfSegments() - 1)[1] + 1 && Coords[2] == GetSegmentCoordinate(GetNumberOfSegments() - 1)[2] + 1)
        //        {
        //            adjacent = true;
        //        }
        //        if (Coords[0] == GetSegmentCoordinate(GetNumberOfSegments() - 1)[0] && Coords[1] == GetSegmentCoordinate(GetNumberOfSegments() - 1)[1] + 1 && Coords[2] == GetSegmentCoordinate(GetNumberOfSegments() - 1)[2])
        //        {
        //            adjacent = true;
        //        }
        //        if (Coords[0] == GetSegmentCoordinate(GetNumberOfSegments() - 1)[0] && Coords[1] == GetSegmentCoordinate(GetNumberOfSegments() - 1)[1] - 1 && Coords[2] == GetSegmentCoordinate(GetNumberOfSegments() - 1)[2])
        //        {
        //            adjacent = true;
        //        }
        //        if (Coords[0] == GetSegmentCoordinate(GetNumberOfSegments() - 1)[0] && Coords[1] == GetSegmentCoordinate(GetNumberOfSegments() - 1)[1] + 1 && Coords[2] == GetSegmentCoordinate(GetNumberOfSegments() - 1)[2] + 1)
        //        {
        //            adjacent = true;
        //        }
        //        if (Coords[0] == GetSegmentCoordinate(GetNumberOfSegments() - 1)[0] && Coords[1] == GetSegmentCoordinate(GetNumberOfSegments() - 1)[1] + 1 && Coords[2] == GetSegmentCoordinate(GetNumberOfSegments() - 1)[2] - 1)
        //        {
        //            adjacent = true;
        //        }
        //
        //        if (Coords[0] == GetSegmentCoordinate(GetNumberOfSegments() - 1)[0] && Coords[1] == GetSegmentCoordinate(GetNumberOfSegments() - 1)[1] - 1 && Coords[2] == GetSegmentCoordinate(GetNumberOfSegments() - 1)[2] - 1)
        //        {
        //            adjacent = true;
        //        }
        //        if (Coords[0] == GetSegmentCoordinate(GetNumberOfSegments() - 1)[0] && Coords[1] == GetSegmentCoordinate(GetNumberOfSegments() - 1)[1] - 1 && Coords[2] == GetSegmentCoordinate(GetNumberOfSegments() - 1)[2] + 1)
        //        {
        //            adjacent = true;
        //        }
        //        if (Coords[0] == GetSegmentCoordinate(GetNumberOfSegments() - 1)[0] && Coords[1] == GetSegmentCoordinate(GetNumberOfSegments() - 1)[1] && Coords[2] == GetSegmentCoordinate(GetNumberOfSegments() - 1)[2] + 1)
        //        {
        //            adjacent = true;
        //        }
        //        if (Coords[0] == GetSegmentCoordinate(GetNumberOfSegments() - 1)[0] && Coords[1] == GetSegmentCoordinate(GetNumberOfSegments() - 1)[1] && Coords[2] == GetSegmentCoordinate(GetNumberOfSegments() - 1)[2] - 1)
        //        {
        //            adjacent = true;
        //        }
        //        if (Coords[0] == GetSegmentCoordinate(GetNumberOfSegments() - 1)[0] + 1 && Coords[1] == GetSegmentCoordinate(GetNumberOfSegments() - 1)[1] && Coords[2] == GetSegmentCoordinate(GetNumberOfSegments() - 1)[2] + 1)
        //        {
        //            adjacent = true;
        //        }
        //        if (Coords[0] == GetSegmentCoordinate(GetNumberOfSegments() - 1)[0] + 1 && Coords[1] == GetSegmentCoordinate(GetNumberOfSegments() - 1)[1] && Coords[2] == GetSegmentCoordinate(GetNumberOfSegments() - 1)[2] - 1)
        //        {
        //            adjacent = true;
        //        }
        //        if (Coords[0] == GetSegmentCoordinate(GetNumberOfSegments() - 1)[0] - 1 && Coords[1] == GetSegmentCoordinate(GetNumberOfSegments() - 1)[1] && Coords[2] == GetSegmentCoordinate(GetNumberOfSegments() - 1)[2] + 1)
        //        {
        //            adjacent = true;
        //        }
        //        if (Coords[0] == GetSegmentCoordinate(GetNumberOfSegments() - 1)[0] - 1 && Coords[1] == GetSegmentCoordinate(GetNumberOfSegments() - 1)[1] && Coords[2] == GetSegmentCoordinate(GetNumberOfSegments() - 1)[2] - 1)
        //        {
        //            adjacent = true;
        //        }
        //
        //
        //        if (!adjacent)
        //        {
        //            cout << "Number of segments: " << GetNumberOfSegments() << "\n";
        //            cout << "Previous segment coordinates: " << GetSegmentCoordinate(GetNumberOfSegments() - 1) << "\n";
        //            cout << "Next segment coordinates: " << Coords << "\n";
        //        }
        //
        //        assert(adjacent);


        mVesselSegmentLocations.push_back(Coords);

        SetLength(0.0);

        for (unsigned i = 1; i < GetNumberOfSegments(); i++)
        {
            double segmentlength = 0;

            if (SPATIAL_DIM == 2)
            {
                segmentlength = pow((pow(double(GetSegmentCoordinate(i)[0] - GetSegmentCoordinate(i - 1)[0]),2) + pow(double(GetSegmentCoordinate(i)[1] - GetSegmentCoordinate(i - 1)[1]),2)),0.5);
            }

            if (SPATIAL_DIM == 3)
            {
                segmentlength = pow((pow(double(GetSegmentCoordinate(i)[0] - GetSegmentCoordinate(i - 1)[0]),2) + pow(double(GetSegmentCoordinate(i)[1] - GetSegmentCoordinate(i - 1)[1]),2) + pow(double(GetSegmentCoordinate(i)[2] - GetSegmentCoordinate(i - 1)[2]),2)),0.5);
            }

            SetLength(GetLength()+segmentlength);
        }
    }
    else
    {
        mVesselSegmentLocations.push_back(Coords);
        SetLength(1);
    }


}

template<unsigned SPATIAL_DIM>
void CaVessel<SPATIAL_DIM>::SetTimeWithLowWallShearStress(double time)
{
    mTimeWithLowWallShearStress = time;
}

template<unsigned SPATIAL_DIM>
void CaVessel<SPATIAL_DIM>::IncrementTimeWithLowWallShearStress(double t)
{
    mTimeWithLowWallShearStress += t;
}

template<unsigned SPATIAL_DIM>
void CaVessel<SPATIAL_DIM>::SetActiveTipCellLocatedAtNode1(bool value)
{
    if (value == true)
    {
        assert(GetNode1()->GetNumberOfAdjoiningVessels() == 1 && GetNode1()->GetAdjoiningVessel(0) == shared());
        assert(mActiveTipCellLocatedAtNode2 == false);
    }
    mActiveTipCellLocatedAtNode1 = value;
}

template<unsigned SPATIAL_DIM>
void CaVessel<SPATIAL_DIM>::SetActiveTipCellLocatedAtNode2(bool value)
{
    if (value == true)
    {
        assert(GetNode2()->GetNumberOfAdjoiningVessels() == 1 && GetNode2()->GetAdjoiningVessel(0) == shared());
        assert(mActiveTipCellLocatedAtNode1 == false);
    }
    mActiveTipCellLocatedAtNode2 = value;
}

template<unsigned SPATIAL_DIM>
void CaVessel<SPATIAL_DIM>::SetRadius(double value)
{
    mRadius = value;
}

template<unsigned SPATIAL_DIM>
void CaVessel<SPATIAL_DIM>::SetPreviousRadius(double value)
{
    mPreviousRadius = value;
}

template<unsigned SPATIAL_DIM>
void CaVessel<SPATIAL_DIM>::SetHaematocritLevel(double value)
{
    mHaematocritLevel = value;
}

template<unsigned SPATIAL_DIM>
void CaVessel<SPATIAL_DIM>::SetFlowVelocity(double value)
{
    mFlowVelocity = value;
}

template<unsigned SPATIAL_DIM>
void CaVessel<SPATIAL_DIM>::SetFlowRate(double value)
{
    mFlowRate = value;
}

template<unsigned SPATIAL_DIM>
void CaVessel<SPATIAL_DIM>::SetImpedance(double value)
{
    mImpedance = value;
}

template<unsigned SPATIAL_DIM>
void CaVessel<SPATIAL_DIM>::SetLength(double value)
{
    mLength = value;
}

template<unsigned SPATIAL_DIM>
void CaVessel<SPATIAL_DIM>::SetWallShearStress(double value)
{
    mWallShearStress = value;
}

template<unsigned SPATIAL_DIM>
void CaVessel<SPATIAL_DIM>::SetViscosity(double value)
{
    mViscosity = value;
}

template<unsigned SPATIAL_DIM>
void CaVessel<SPATIAL_DIM>::SetMechanicalStimulus(double value)
{
    mMechanicalStimulus = value;
}

template<unsigned SPATIAL_DIM>
void CaVessel<SPATIAL_DIM>::SetMetabolicStimulus(double value)
{
    mMetabolicStimulus = value;
}

template<unsigned SPATIAL_DIM>
void CaVessel<SPATIAL_DIM>::SetShrinkingStimulus(double value)
{
    mShrinkingStimulus = value;
}

template<unsigned SPATIAL_DIM>
void CaVessel<SPATIAL_DIM>::SetDownstreamConductedStimulus(double value)
{
    mDownstreamConductedStimulus = value;
}

template<unsigned SPATIAL_DIM>
void CaVessel<SPATIAL_DIM>::SetUpstreamConductedStimulus(double value)
{
    mUpstreamConductedStimulus = value;
}

template<unsigned SPATIAL_DIM>
bool CaVessel<SPATIAL_DIM>::IsAttachedToNode(boost::shared_ptr<CaVascularNetworkNode<SPATIAL_DIM> > node)
{
    return (node == pNode1 || node == pNode2);
}

template<unsigned SPATIAL_DIM>
void CaVessel<SPATIAL_DIM>::SetIsPartOfNeovasculature(bool value)
{
    mIsPartOfNeovasculature = value;
}

template<unsigned SPATIAL_DIM>
void CaVessel<SPATIAL_DIM>::CanExtend(bool value)
{
    mCanExtend = value;
}



// Explicit instantiation

template class CaVessel<2>;
template class CaVessel<3>;


