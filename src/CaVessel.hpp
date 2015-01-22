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

#ifndef CAVESSEL_HPP_
#define CAVESSEL_HPP_

#include <vector>
#include <iostream>
#include <map>
#include <boost/enable_shared_from_this.hpp>
#include "SmartPointers.hpp"
#include "Exception.hpp"
#include "ChastePoint.hpp"
#include "IntraVascularChemicalCollection.hpp"
#include "Concentration.hpp"
#include "CaVascularNetworkNode.hpp"

template<unsigned DIM>
class CaVessel : public boost::enable_shared_from_this<CaVessel<DIM> >
{

private:

	/**
	 *  Node at one end of the Vessel.
	 */
	boost::shared_ptr<CaVascularNetworkNode<DIM> > pNode1;

	/**
	 *  Node at other end of the Vessel.
	 */
	boost::shared_ptr<CaVascularNetworkNode<DIM> > pNode2;

	/**
	 *   Vessel data of type double, with units
	 */
	std::map<std::string, std::pair<double, std::string> > mDoubleData;

	/**
	 *   Vessel data of type bool
	 */
	std::map<std::string, bool> mBooleanData;

	/**
	 *  Collection of VesselSegments. A vessel segment is a spatially discrete location
	 *  occupied by the Vessel.
	 *
	 *  todo Need to think about how we introduce the following with cells instead of segments. For now use a vector of locations.
	 */
	std::vector<ChastePoint<DIM> > mVesselSegmentLocations;

	/**
	 * Collection of chemicals which the Vessel contains.
	 */
	IntraVascularChemicalCollection mChemicalCollection;

public:

	/**
       Constructor.

       When a CaVessel object is instantiated the associated CaVascularNetworkNodes have no location.
       These must be prescribed before a CaVessel is added to a CaVesselNetwork. The collection of
       segment coordinates is also empty. When segement coordinates are added to the collection of
       vessel segement coordinates collection, they must be added in order. Each coordinate must be
       located at a lattice site which neighbours the location of the previous segment coordinate
       added to the collection. A CaVessel initially contains no IntraVascularChemicals.
	 */
	CaVessel();

	/**
       Destructor.
	 */
	~CaVessel();

	/**
       Copies the mechanical property values and chemical concentrations from the prescribed vessel
       to this vessel. This method is useful for manipulating vessels within a CaVesselNetwork object.
       For example, this method is useful for when an existing vessel divides to form two new vessels.
       None of the spatial information about the prescribed vessel are copied (i.e. segment coordinates).
	 */
	void CopyMechanicalPropertyValuesAndChemicalConcentrations(boost::shared_ptr<CaVessel<DIM> > anotherVessel);

	/**
       @return boost::shared_ptr<CaVessel<DIM> >
	 */
	boost::shared_ptr<CaVessel<DIM> > shared();

	/**
       @return shared pointer to node1.
	 */
	boost::shared_ptr<CaVascularNetworkNode<DIM> > GetNode1();

	/**
       @return shared pointer to node2.
	 */
	boost::shared_ptr<CaVascularNetworkNode<DIM> > GetNode2();

	/**
       @return mVesselSegmentLocations
	 */
	std::vector<ChastePoint<DIM> > GetSegmentCoordinates();

	/**
       @return mVesselSegmentLocations[i]
	 */
	ChastePoint<DIM> GetSegmentCoordinate(unsigned i);

	/**
       @return mVesselSegmentLocations.size()
	 */
	unsigned GetNumberOfSegments();

	/**
	 * Returns type double vessel data value.
	 */
	double GetDoubleDataValue(const std::string& variableName);

	/**
	 * Returns type double vessel data units.
	 */
	const std::string& GetDoubleDataUnits(const std::string& variableName) const;

	/**
	 * Returns type boolean vessel data.
	 */
	bool GetBooleanData(const std::string& variableName);

	/**
	 * Assigns type double vessel data.
	 */
	void SetDoubleData(const std::string& variableName, double data, const std::string& unit = "None");

	/**
	 * Assigns type boolean vessel data.
	 */
	void SetBooleanData(const std::string& variableName, bool data);

	/**
	 *  @return whether vessel is attached to input node.
	 */
	bool IsInputVessel();

	/*
	 * todo A variant of the following should be added when cells (rather than vessel segments) are associated with vessels.
	 */
	//    /**
	//        @return mHaematocrit
	//     */
	//    boost::shared_ptr<VesselSegment<DIM> > GetVesselSegment(int i);
	//
	//    /**
	//         Returns the vessel segment at the specified index
	//     */
	//    boost::shared_ptr<VesselSegment<DIM> > GetVesselSegment(ChastePoint<DIM> loc);

	/**
       @return mChemicalCollection
	 */
	IntraVascularChemicalCollection& rGetCollectionOfIntraVascularChemicals();

	/**
       @return the number of IntraVascularChemicals contained within the vessel.
	 */
	unsigned GetNumberOfIntraVascularChemicals();

	/**
       @return the concentration of the prescribed IntraVascularChemical.
	 */
	double GetIntraVascularChemicalConcentration(string chemicalname);

	/**
       @return the tortuosity of a vessel. The tortuosity is calculated as the arc-chord ratio. This is the ratio of the vessel length to the distance between the two ends of the vessel.
	 */
	double GetTortuosity();

	/**
       Assigns the concentration of the prescribed IntraVascularChemical.

       @param chemicalname
       @param concentration
	 */
	void SetIntraVascularChemicalConcentration(string chemicalname, Concentration concentration);

	/**
       Assigns the prescribed node to pNode1;

       @param node node to be assigned to pNode1
	 */
	void SetNode1(boost::shared_ptr<CaVascularNetworkNode<DIM> > node);

	/**
       Assigns the prescribed location to pNode1 of the vessel.

       todo Within this method the vessel is also added to the node as an adjoining vessel. Should
       look in to a way of improving how adjoining vessels are added to nodes and nodes to vessels,
       whilst maintaining a consistent state in the system.

       @param location new location of pNode1
	 */
	void SetNode1Location(ChastePoint<DIM> location);

	/**
       Assigns the prescribed node to pNode2

       @param node node to be assigned to pNode2
	 */
	void SetNode2(boost::shared_ptr<CaVascularNetworkNode<DIM> > node);

	/**
       Assigns the prescribed location to pNode2 of the vessel.

       todo Within this method the vessel is also added to the node as an adjoining vessel. Should
       look in to a way of improving how adjoining vessels are added to nodes and nodes to vessels,
       whilst maintaining a consistent state in the system.

       @param location new location of pNode2
	 */
	void SetNode2Location(ChastePoint<DIM> location);

	/**
       Assigns a new vessel segment location to the collection of vessel segments. Vessel segments
       must be added in order and each new vessel segment must be located at a lattice site which
       neighbours the last vessel segment added to the selection.

       This method also re-calculates the length of the vessel, given that it has had
       another segment added to it.

       todo We could also potentially also automatically update the location of CaVascularNetworkNodes
       in this method. Would avoid manually updating node locations after all vessel segments have
       been added.

       @param location location of next vessel segment
	 */
	void SetNextVesselSegmentCoordinate(ChastePoint<DIM> location);

	/*
	 * todo Potentially add the following when cells are explicitly associated with vessels.
	 */
	//    /**
	//            Add new VesselSegment to vessel
	//     */
	//    void AddVesselSegment(ChastePoint<DIM> Coords);

	/**
       Checks whether the prescribed node is attached to this vessel.

       @return true if vessel is attached to prescribed node
	 */
	bool IsAttachedToNode(boost::shared_ptr<CaVascularNetworkNode<DIM> > node);

};

#endif /* CAVESSEL_HPP_ */
