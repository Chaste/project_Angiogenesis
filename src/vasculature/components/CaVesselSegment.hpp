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

#ifndef CAVESSELSEGMENT_HPP_
#define CAVESSELSEGMENT_HPP_

#include <boost/enable_shared_from_this.hpp>
#include <math.h>
#include "VascularNode.hpp"
#include "CaVessel.hpp"
#include "VasculatureData.hpp"
#include "SmartPointers.hpp"
#include "ChastePoint.hpp"

template<unsigned DIM>
class CaVessel;

template<unsigned DIM>
class VascularNode;

/*
 * Vessel segments are straight sub-units of vessels, defined by the positions of
 * their end nodes. Nodes cannot be created by the vessel segment class, they are
 * instead managed by the VascularNetwork class. Segments must always have two nodes.
 */

template<unsigned DIM>
class CaVesselSegment : public boost::enable_shared_from_this<CaVesselSegment<DIM> >
{

	///\todo this breaks encapsulation, but ensures that the segment-vessel connectivity is
	// kept up-to-date. Look at methods for limiting friend class access to methods.
	friend class CaVessel<DIM>;

private:

    /**
     * Container for segment end nodes
     */
    std::pair<boost::shared_ptr<VascularNode<DIM> >, boost::shared_ptr<VascularNode<DIM> > > mNodes;

    /**
     * Container for non-spatial segment data.
     */
	boost::shared_ptr<VasculatureData> mpDataContainer;

    /**
     * Id tag, can be useful for storing segment-vessel relationships in the VesselNetwork class.
     */
    unsigned mId;

    /**
     * Label tag, can be useful for identifying input and output segments.
     */
    std::string mLabel;

	/**
	 *   Weak pointer to the vessel owning this segment
	 */
	boost::weak_ptr<CaVessel<DIM> > mVessel;

private:

    /*
     * Constructor - This is private as instances of this class must be created with a corresponding shared pointer. This is
     * implemented using the static Create method.
     */
    CaVesselSegment(boost::shared_ptr<VascularNode<DIM> > pNode1, boost::shared_ptr<VascularNode<DIM> > pNode2);

public:

    /*
     * Construct a new instance of the class and return a shared pointer to it. Also manage the association of segments to nodes by
     * passing weak pointers to the nodes.
     */
    static boost::shared_ptr<CaVesselSegment<DIM> > Create(boost::shared_ptr<VascularNode<DIM> > pNode1, boost::shared_ptr<VascularNode<DIM> > pNode2);

    /*
     * Destructor
     */
    ~CaVesselSegment();

	/**
	 *  Return a pointer to the segment's non-spatial data container.
	 *
	 *  @return mDataContainer
	 */
	boost::shared_ptr<VasculatureData> GetDataContainer();

	/**
	 *  Return the Id
	 *
	 *  @return mId
	 */
	unsigned GetId();

	/**
	 *  Return a boost::shared_ptr to the vessel.
	 *
	 * @return mAdjoiningVessels[i]
	 */
	boost::shared_ptr<CaVessel<DIM> > GetVessel();

	/**
	 *  Return the Label
	 *
	 *  @return mLabel
	 */
	const std::string& rGetLabel();

	/**
	 *  Return the length
	 */
    double GetLength();

	/**
	 *  Return a point mid-way along the vessel
	 */
    ChastePoint<DIM> GetMidPoint();

    /**
       Return a the segment nodes as a pair
     */
    std::pair<boost::shared_ptr<VascularNode<DIM> >, boost::shared_ptr<VascularNode<DIM> > > GetNodes();

    /**
       Return a pointer to the node specified by the index
     */
    boost::shared_ptr<VascularNode<DIM> > GetNodes(unsigned index);

    /**
       Replace the node at the specified index with the passed in node.
     */
    void ReplaceNode(unsigned old_node_index, boost::shared_ptr<VascularNode<DIM> >  pNewNode);

	/**
	 *  Over-write the segment's non-spatial DataContainer
	 *
	 *  This can be useful when copying data from an existing segment.
	 */
	void SetDataContainer(boost::shared_ptr<VasculatureData> pDataContainer);

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

private:

    /**
       Return a boost::shared_ptr to this object

       @return boost::shared_ptr<VesselSegment<DIM> >
     */
	boost::shared_ptr<CaVesselSegment<DIM> > Shared();

	/**
       Adds an adjoining Vessel to the segment.
	 */
	void AddVessel(boost::shared_ptr<CaVessel<DIM> > pVessel);

	/**
       Removes an adjoining vessel from the segment.
	 */
	void RemoveVessel();
};

#endif /* CAVASCULARNETWORK_HPP_ */
