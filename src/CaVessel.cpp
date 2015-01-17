/*
 * CaVessel.cpp
 *
 *  Created on: 13 Jan 2015
 *      Author: connor
 */

#include "CaVessel.hpp"

template<unsigned DIM>
CaVessel<DIM>::CaVessel() : pNode1(new CaVascularNetworkNode<DIM>()),pNode2(new CaVascularNetworkNode<DIM>()),mActiveTipCellLocatedAtNode1(false),mActiveTipCellLocatedAtNode2(false),mVesselSegmentLocations(),mTimeWithLowWallShearStress(0.0),mRadius(1.0*pow(10.0,-5)),mPreviousRadius(1.0),mHaematocritLevel(0.45),mFlowVelocity(0.0),mFlowRate(0.0),mImpedance(0.0),mLength(0.0),mWallShearStress(0.0),mViscosity(0.0),mMechanicalStimulus(0.0),mMetabolicStimulus(0.0),mShrinkingStimulus(0.0),mDownstreamConductedStimulus(0.0),mUpstreamConductedStimulus(0.0),mChemicalCollection(),mIsPartOfNeovasculature(true), mCanExtend(false)
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
void CaVessel<DIM>::CopyMechanicalPropertyValuesAndChemicalConcentrations(boost::shared_ptr<CaVessel<DIM> > anotherVessel)
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

template<unsigned DIM>
double CaVessel<DIM>::GetRadius()
{
    return mRadius;
}

template<unsigned DIM>
double CaVessel<DIM>::GetTortuosity()
{
    double distanceBetweenEndsOfVessel;

    if (DIM > 2)
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
double CaVessel<DIM>::GetTimeWithLowWallShearStress() const
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
    assert(i < mVesselSegmentLocations.size());
    return mVesselSegmentLocations[i];
}

template<unsigned DIM>
unsigned CaVessel<DIM>::GetNumberOfSegments()
{
    return mVesselSegmentLocations.size();
}


template<unsigned DIM>
bool CaVessel<DIM>::HasActiveTipCell() const
{
    return (ActiveTipCellLocatedAtNode1() || ActiveTipCellLocatedAtNode2());
}

template<unsigned DIM>
bool CaVessel<DIM>::ActiveTipCellLocatedAtNode1() const
{
    return mActiveTipCellLocatedAtNode1;
}

template<unsigned DIM>
bool CaVessel<DIM>::ActiveTipCellLocatedAtNode2() const
{
    return mActiveTipCellLocatedAtNode2;
}

template<unsigned DIM>
IntraVascularChemicalCollection& CaVessel<DIM>::GetCollectionOfIntraVascularChemicals()
{
    return mChemicalCollection;
}

template<unsigned DIM>
unsigned CaVessel<DIM>::GetNumberOfIntraVascularChemicals()
{
    return mChemicalCollection.GetIntraVascularChemicalCollection().size();
}

template<unsigned DIM>
double CaVessel<DIM>::GetIntraVascularChemicalConcentration(string chemicalname)
{
    return mChemicalCollection.GetIntraVascularChemicalConcentration(chemicalname);
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
void CaVessel<DIM>::SetIntraVascularChemicalConcentration(string chemicalname, Concentration concentration)
{
    mChemicalCollection.SetIntraVascularChemicalConcentration(chemicalname,concentration);
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
 * todo This method needs some adjustment ... dealing with coordinates that may not be integers, for example. And what if there is no z-coordinate
 */
template<unsigned DIM>
void CaVessel<DIM>::SetNextVesselSegmentCoordinate(ChastePoint<DIM> Coords)
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

            if (DIM == 2)
            {
                segmentlength = pow((pow(double(GetSegmentCoordinate(i)[0] - GetSegmentCoordinate(i - 1)[0]),2) + pow(double(GetSegmentCoordinate(i)[1] - GetSegmentCoordinate(i - 1)[1]),2)),0.5);
            }

            if (DIM == 3)
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

template<unsigned DIM>
void CaVessel<DIM>::SetTimeWithLowWallShearStress(double time)
{
    mTimeWithLowWallShearStress = time;
}

template<unsigned DIM>
void CaVessel<DIM>::IncrementTimeWithLowWallShearStress(double t)
{
    mTimeWithLowWallShearStress += t;
}

template<unsigned DIM>
void CaVessel<DIM>::SetActiveTipCellLocatedAtNode1(bool value)
{
    if (value == true)
    {
        assert(GetNode1()->GetNumberOfAdjoiningVessels() == 1 && GetNode1()->GetAdjoiningVessel(0) == shared());
        assert(mActiveTipCellLocatedAtNode2 == false);
    }
    mActiveTipCellLocatedAtNode1 = value;
}

template<unsigned DIM>
void CaVessel<DIM>::SetActiveTipCellLocatedAtNode2(bool value)
{
    if (value == true)
    {
        assert(GetNode2()->GetNumberOfAdjoiningVessels() == 1 && GetNode2()->GetAdjoiningVessel(0) == shared());
        assert(mActiveTipCellLocatedAtNode1 == false);
    }
    mActiveTipCellLocatedAtNode2 = value;
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


