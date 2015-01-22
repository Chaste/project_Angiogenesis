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
	  mDoubleData(),
	  mBooleanData(),
	  mVesselSegmentLocations(),
	  mChemicalCollection()
{
}

template<unsigned DIM>
CaVessel<DIM>::~CaVessel()
{
}

//template<unsigned DIM>
//bool CaVessel<DIM>::IsInputVessel()
//{
//    return (pNode1->IsInputNode() || pNode2->IsInputNode());
//}

template<unsigned DIM>
void CaVessel<DIM>::CopyMechanicalPropertyValuesAndChemicalConcentrations(boost::shared_ptr<CaVessel<DIM> > another_vessel)
{
    // todo need to copy property maps over

    for (unsigned i = 0; i < another_vessel->GetNumberOfIntraVascularChemicals(); i++)
    {
        mChemicalCollection.AddIntraVascularChemical(another_vessel->rGetCollectionOfIntraVascularChemicals().GetIntraVascularChemicalCollection()[i].GetChemicalName(),
                Concentration(another_vessel->rGetCollectionOfIntraVascularChemicals().GetIntraVascularChemicalCollection()[i].GetConcentration(),
                another_vessel->rGetCollectionOfIntraVascularChemicals().GetIntraVascularChemicalCollection()[i].GetUnits()),
                another_vessel->rGetCollectionOfIntraVascularChemicals().GetIntraVascularChemicalCollection()[i].GetPermeabilityOfVesselWallToChemical());
    }
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

    // todo need to change the following to have length returned from property list

    double length = 100;

    return (length/distance_between_ends_of_vessel);
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
double CaVessel<DIM>::GetDoubleDataValue(const std::string& variableName)
{
	std::map<std::string, std::pair<double, std::string> >::const_iterator it = mDoubleData.find(variableName);
	if (it == mDoubleData.end())
	{
		EXCEPTION("No double valued property, '" << variableName << "', in property register.");
	}
	return(it->second.first);
}

template<unsigned DIM>
const std::string& CaVessel<DIM>::GetDoubleDataUnits(const std::string& variableName) const
{

	std::map<std::string, std::pair<double, std::string> >::const_iterator it = mDoubleData.find(variableName);
	if (it == mDoubleData.end())
	{
		EXCEPTION("No double valued property, '" << variableName << "', in property register.");
	}
	return(it->second.second);

}

template<unsigned DIM>
bool CaVessel<DIM>::GetBooleanData(const std::string& variableName)
{
	std::map<std::string, bool >::const_iterator it = mBooleanData.find(variableName);
	if (it == mBooleanData.end())
	{
		EXCEPTION("No boolean valued property, '" << variableName << "', in property register.");
	}
	return(it->second);
}

template<unsigned DIM>
void CaVessel<DIM>::SetDoubleData(const std::string& variableName, double data, const std::string& unit)
{
	mDoubleData[variableName] = std::pair<double, std::string> (data, unit);
}

template<unsigned DIM>
void CaVessel<DIM>::SetBooleanData(const std::string& variableName, bool data)
{
	mBooleanData[variableName] = data;
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
//    if (GetNumberOfSegments() > 0)
//    {
//        mVesselSegmentLocations.push_back(point);
//        // todo set length here?
//
//        for (unsigned i = 1; i < GetNumberOfSegments(); i++)
//        {
//            double segment_length = 0;
//
//            if (DIM == 2)
//            {
//                segment_length = pow((pow(double(GetSegmentCoordinate(i)[0]-GetSegmentCoordinate(i-1)[0]), 2) +
//                		pow(double(GetSegmentCoordinate(i)[1]-GetSegmentCoordinate(i-1)[1]), 2)), 0.5);
//            }
//
//            if (DIM == 3)
//            {
//                segment_length = pow((pow(double(GetSegmentCoordinate(i)[0]-GetSegmentCoordinate(i-1)[0]), 2) +
//                		pow(double(GetSegmentCoordinate(i)[1]-GetSegmentCoordinate(i-1)[1]),2) +
//                		pow(double(GetSegmentCoordinate(i)[2] - GetSegmentCoordinate(i - 1)[2]), 2)), 0.5);
//            }
//            // todo set length here?
//        }
//    }
//    else
//    {
//        mVesselSegmentLocations.push_back(point);
//        // todo set length here?
//    }
}

template<unsigned DIM>
bool CaVessel<DIM>::IsAttachedToNode(boost::shared_ptr<CaVascularNetworkNode<DIM> > node)
{
    return (node == pNode1 || node == pNode2);
}

// Explicit instantiation
template class CaVessel<2>;
template class CaVessel<3>;
