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
#include <algorithm>
#include <boost/weak_ptr.hpp>
#include "ChastePoint.hpp"
#include "Exception.hpp"
#include "SmartPointers.hpp"
#include "Cell.hpp"
#include "CaBasedCellPopulation.hpp"
#include "VascularNetworkData.hpp"
#include "CaVesselSegment.hpp"

template<unsigned DIM>
class CaVesselSegment;

/*
 * Nodes are point locations along a vessel. They are useful for describing the positions of
 * vessel segments and in calculating flow in networks. They can be linked to the position of a single Cell object
 * which can in turn correspond to a single biological cell or a group of biological cells. Cell objects should update
 * nodal positions while nodes can pass flow and stimulus information to their corresponding Cells.
 */

template<unsigned DIM>
class CaVascularNetworkNode : public boost::enable_shared_from_this<CaVascularNetworkNode<DIM> >
{
	///\todo this breaks encapsulation, but ensures that the node-segment connectivity is
	// kept up-to-date. Look at methods for limiting friend class access to methods.
	friend class CaVesselSegment<DIM>;

private:

	/**
	 *   Location of a node. Used if a CellPtr has not been assigned.
	 */
	ChastePoint<DIM> mLocation;

	/**
	 *   Pointer to an associated Cell.
	 */
	CellPtr mpCell;

	/**
	 *   Pointer to an associated Cell Population.
	 */
	CaBasedCellPopulation<DIM>* mpCellPopulation;

    /**
     * Container for non-spatial node data.
     */
	boost::shared_ptr<VascularNetworkData> mpDataContainer;

    /**
     * Id tag, can be useful for storing segment-node relationships in the VesselNetwork class.
     */
    unsigned mId;

    /**
     * Label tag, can be useful for identifying input and output nodes.
     */
    std::string mLabel;

	/**
	 *   Collection of weak pointers to Vessel Segments connected to this node.
	 */
	std::vector<boost::weak_ptr<CaVesselSegment<DIM> > > mVesselSegments;

public:

	/*
	 * Constructor
	 *
	 * A node must be specified with a location.
	 */
	CaVascularNetworkNode(ChastePoint<DIM> location);

    /*
     * Construct a new instance of the class and return a shared pointer to it.
     */
    static boost::shared_ptr<CaVascularNetworkNode<DIM> > Create(ChastePoint<DIM> location);

	/*
	 * Destructor
	 */
	~CaVascularNetworkNode();

	/**
	 *  Return a pointer to the associated Cell.
	 *
	 *  @return mpCell
	 */
	CellPtr GetCell();

	/**
	 *  Return the Id
	 *
	 *  @return mId
	 */
	unsigned GetId();

	/**
	 *  Return the Label
	 *
	 *  @return mLabel
	 */
	const std::string& rGetLabel();

	/**
	 *  Return the location of the node or, if there is one, the associated Cell.
	 *
	 *  @return mLocation
	 */
	ChastePoint<DIM> GetLocation();

	/**
	 *  Return a pointer to the node's non-spatial data container.
	 *
	 *  @return mDataContainer
	 */
	boost::shared_ptr<VascularNetworkData> GetDataContainer();

	/**
	 *  Return the number of attached segments
	 */
	unsigned GetNumberOfSegments();

	/**
	 *  Return a boost::shared_ptr to VesselSegment i.
	 *
	 * @return mAdjoiningVessels[i]
	 */
	boost::shared_ptr<CaVesselSegment<DIM> > GetVesselSegments(unsigned index);

	std::vector<boost::shared_ptr<CaVesselSegment<DIM> > >GetVesselSegments();

	/**
	 *  Return true if there is an associated Cell.
	 *
	 *  @return bool
	 */
	bool HasCell();

	/**
	 *  Assign a Cell to the node. Throw an Exception if the Cell is not a member of the
	 *  Cell Population or if a Cell Population has not been set. Overwrite any existing Cell.
	 */
	void SetCell(CellPtr pCell);

	/**
	 *  Assign a Cell Population to the node. If an existing Cell is not a member of the population
	 *  remove it.
	 */
	void SetCellPopulation(CaBasedCellPopulation<DIM>* pCellPopulation);

	/**
	 *  Over-write the node's non-spatial DataContainer
	 *
	 *  This can be useful when copying data from an existing node.
	 */
	void SetDataContainer(boost::shared_ptr<VascularNetworkData> pDataContainer);

	/**
	 *  Assign the Id
	 *
	 */
	void SetId(unsigned id);

	/**
	 *  Assign the Label
	 *
	 */
	void SetLabel(const std::string& label);

	/**
	 *  Set the location of the node. This breaks any links with an assigned Cell, so if there is an
	 *  assigned Cell remove it.
	 */
	void SetLocation(ChastePoint<DIM> location);

	/**
	 *  Remove the assigned Cell
	 */
	void RemoveCell();

private:

	/**
       Adds an adjoining Vessel segment the node.
	 */
	void AddSegment(boost::shared_ptr<CaVesselSegment<DIM> > pVesselSegment);

	/**
       Removes an adjoining vessel segment from to node.
	 */
	void RemoveSegment(boost::shared_ptr<CaVesselSegment<DIM> > pVesselSegment);

};

#endif /* CAVASCULARNETWORKNODE_HPP_ */
