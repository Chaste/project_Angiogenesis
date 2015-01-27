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
#include "CaVascularNetworkNode.hpp"
#include "CaVesselSegment.hpp"

template<unsigned DIM>
class CaVesselSegment;

template<unsigned DIM>
class CaVascularNetworkNode;

template<unsigned DIM>
class CaVessel : public boost::enable_shared_from_this<CaVessel<DIM> >
{
private:

	/**
	 *  Vessel segments
	 */
	std::vector<boost::shared_ptr<CaVesselSegment<DIM> > > mSegments;

    /**
     * Container for non-spatial vessel data.
     */
	boost::shared_ptr<VascularNetworkData> mpDataContainer;

    /**
     * Id tag, can be useful for storing segment-vessel relationships in the VesselNetwork class.
     */
    unsigned mId;

    /**
     * Label tag, can be useful for identifying input and output vessels.
     */
    std::string mLabel;

private:

	/**
       Constructor.

       The vessel should always have at least one segment.
	 */
	CaVessel(boost::shared_ptr<CaVesselSegment<DIM> > pSegment);

	/**
       Alternate Constructor.

       The vessel should always have at least one segment. This is useful for initializing with many segments at once.
       ///\todo initializing with a list of nodes would also be handy.
	 */
	CaVessel(std::vector<boost::shared_ptr<CaVesselSegment<DIM> > > segments);

public:

    /*
     * Construct a new instance of the class and return a shared pointer to it. Also manage the association of segments to nodes by
     * passing weak pointers to the nodes.
     */
    static boost::shared_ptr<CaVessel<DIM> > Create(boost::shared_ptr<CaVesselSegment<DIM> > pSegment);

    /*
     * Construct a new instance of the class and return a shared pointer to it. Also manage the association of segments to nodes by
     * passing weak pointers to the nodes.
     */
    static boost::shared_ptr<CaVessel<DIM> > Create(std::vector<boost::shared_ptr<CaVesselSegment<DIM> > > segments);

	/**
       Destructor.
	 */
	~CaVessel();

	/**
       Add a single segment to either end of the vessel
	 */
	void AddSegments(boost::shared_ptr<CaVesselSegment<DIM> > segment);

	/**
       Add a collection of segments to either end of the vessel
	 */
	void AddSegments(std::vector<boost::shared_ptr<CaVesselSegment<DIM> > > segments);

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
	 *  Return a pointer to the vessel's non-spatial data container.
	 *
	 *  @return mDataContainer
	 */
	boost::shared_ptr<VascularNetworkData> GetDataContainer();

	/**
	 *  Over-write the vessel's non-spatial DataContainer
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
       @return boost::shared_ptr<CaVessel<DIM> >
	 */
	boost::shared_ptr<CaVessel<DIM> > Shared();

	/**
       @return shared pointer to the first node of the first segment
	 */
	boost::shared_ptr<CaVascularNetworkNode<DIM> > GetStartNode();

	/**
       @return shared pointer to the second node of the last segment
	 */
	boost::shared_ptr<CaVascularNetworkNode<DIM> > GetEndNode();

	/**
       @return mVesselSegmentLocations
	 */
	std::vector<boost::shared_ptr<CaVesselSegment<DIM> > > GetSegments();

	/**
       @return mVesselSegmentLocations[i]
	 */
	boost::shared_ptr<CaVesselSegment<DIM> > GetSegments(unsigned i);

	/**
       @return mVesselSegmentLocations.size()
	 */
	unsigned GetNumberOfSegments();
};

#endif /* CAVESSEL_HPP_ */
