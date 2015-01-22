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

#ifndef CAVASCULARNETWORKNODE_HPP_
#define CAVASCULARNETWORKNODE_HPP_

#include <vector>
#include <map>
#include <string>
#include "ChastePoint.hpp"
#include "Exception.hpp"
#include "SmartPointers.hpp"

// forward declaration of vessel class - vascular network node referenced in vessel class
template <unsigned DIM>
class CaVessel;

template<unsigned DIM>
class CaVascularNetworkNode : public boost::enable_shared_from_this<CaVascularNetworkNode<DIM> >
{
private:

	/**
	 *   Location of node.
	 */
	ChastePoint<DIM> mLocation;

	/**
	 *   Node data of type double, with units
	 */
	std::map<std::string, std::pair<double, std::string> > mDoubleData;

	/**
	 *   Node data of type bool
	 */
	std::map<std::string, bool> mBooleanData;

	/**
	 *   Collection of CaVessel objects which are adjoint to this node.
	 */
	std::vector<boost::weak_ptr<CaVessel<DIM> > > mAdjoiningVessels;

public:

	/*
	 * Constructor
	 *
	 * Upon instantiation the node has no prescribed location or adjoining vessels.
	 *
	 */
	CaVascularNetworkNode();

	/*
	 * Destructor
	 */
	~CaVascularNetworkNode();

	/**
	 *  Returns a boost::shared_ptr to this CaVascularNetworkNode.
	 *
	 *  @return boost::shared_ptr<CaVascularNetworkNode<DIM,T> >
	 */
	boost::shared_ptr<CaVascularNetworkNode<DIM> > shared();

	/**
	 *  Returns the location of the node.
	 *
	 *  @return mLocation
	 */
	ChastePoint<DIM> GetLocation();

	/**
	 * Returns type double vessel node data value.
	 */
	double GetDoubleDataValue(const std::string& variableName);

	/**
	 * Returns type double vessel node data units.
	 */
	const std::string& GetDoubleDataUnits(const std::string& variableName) const;

	/**
	 * Returns type boolean vessel node data.
	 */
	bool GetBooleanData(const std::string& variableName);

	/**
	 *  Returns the number of vessels which are adjoint to the node.
	 *
	 *  @return mAdjoiningVessels.size()
	 */
	unsigned GetNumberOfAdjoiningVessels();

	/**
	 *  Returns a boost::shared_ptr to Vessel i which is adjoint to this node.
	 *
	 * @return mAdjoiningVessels[i]
	 */
	boost::shared_ptr<CaVessel<DIM> > GetAdjoiningVessel(unsigned i);

	/**
	 * Assigns type double vessel node data.
	 */
	void SetDoubleData(const std::string& variableName, double data, const std::string& unit = "None");

	/**
	 * Assigns type boolean vessel node data.
	 */
	void SetBooleanData(const std::string& variableName, bool data);

	/**
	 * Assigns the prescribed location to the node.
	 *
	 *  @param loc new location of node
	 */
	void SetLocation(ChastePoint<DIM> location);

	/**
       Adds an adjoining Vessel to the node.
       The prescribed Vessel may be attached to the same node twice if it loops around on itself.
       In this case the node must be both node1 and node2 of the Vessel.

       todo: at some point should possibly check that the vessel being adjoined to node has not got
       an actively migrating tip located at this node.

       @param vessel vessel to be added to node
	 */
	void AddAdjoiningVessel(boost::shared_ptr<CaVessel<DIM> > vessel);

	/**
       Removes an adjoining vessel from to node.

       @param vessel vessel to be removed from node
	 */
	void RemoveAdjoiningVessel(boost::shared_ptr<CaVessel<DIM> > vessel);

	/**
       Checks whether the prescribed vessel is attached to this node.

       @param vessel vessel which may or may not be attached to node
	 */
	bool IsAttachedToVessel(boost::shared_ptr<CaVessel<DIM> > vessel);

};

#endif /* CAVASCULARNETWORKNODE_HPP_ */
