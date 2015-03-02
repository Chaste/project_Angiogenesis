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
#include "SmartVasculaturePointers.hpp"
#include "Exception.hpp"
#include "ChastePoint.hpp"
#include "VascularNode.hpp"
#include "CaVesselSegment.hpp"
#include "VasculatureData.hpp"

struct SegmentLocation
{
    enum Value
    {
        Start,
        End
    };
};

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
	VasculatureData mDataContainer;

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
	 */
	CaVessel(std::vector<boost::shared_ptr<CaVesselSegment<DIM> > > segments);
	
	/**
       Alternate Constructor.

       Initialize with a vector of nodes. The nodes are joined by segments in order. The ends are not closed.
	 */
	CaVessel(std::vector<boost::shared_ptr<VascularNode<DIM> > > nodes);

public:

    /*
     * Construct a new instance of the class and return a shared pointer to it.
     */
    static boost::shared_ptr<CaVessel<DIM> > Create(boost::shared_ptr<CaVesselSegment<DIM> > pSegment);

    /*
     * Construct a new instance of the class and return a shared pointer to it.
     */
    static boost::shared_ptr<CaVessel<DIM> > Create(std::vector<boost::shared_ptr<CaVesselSegment<DIM> > > segments);
    
    /*
     * Construct a new instance of the class and return a shared pointer to it.
     */
    static boost::shared_ptr<CaVessel<DIM> > Create(std::vector<boost::shared_ptr<VascularNode<DIM> > > nodes);

	/**
       Destructor.
	 */
	~CaVessel();

	/**
       Add a single segment to either end of the vessel
	 */
	void AddSegment(boost::shared_ptr<CaVesselSegment<DIM> > pSegment);

	/**
       Add a collection of segments to either end of the vessel
	 */
	void AddSegments(std::vector<boost::shared_ptr<CaVesselSegment<DIM> > > pSegments);
	
    /**
     *  Return the vessel data for the input key. An attempt is made
     *  to cast to type T.
     *  @return T data
     */
    template<typename T> T GetData(const std::string& rKey);

	/**
	 *  Return a const reference to the vessel's non-spatial data container.
	 *
	 *  @return mDataContainer
	 */
	const VasculatureData& rGetDataContainer() const;

    /**
     *  Return a vector of data keys for the vessel. Input true if
     *  the corresponding value should be castable to double.
     *
     *  @return std::vector<std::string>
     */
    std::vector<std::string> GetDataKeys(bool castable_to_double = false) const;

	/**
       @return shared pointer to the second node of the last segment
	 */
	boost::shared_ptr<VascularNode<DIM> > GetEndNode();

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
	 *  Return the length
	 *
	 *  @return double
	 */
	double GetLength();

	/**
	 *  Return the vessel's nodes
	 *
	 *  @return mLabel
	 */
	std::set<boost::shared_ptr<VascularNode<DIM> > > GetNodes();
	
	/**
        Return the nodes in the vessel in the form of a vector
	 */
	std::vector<boost::shared_ptr<VascularNode<DIM> > > GetVectorOfNodes();

	/**
	 * Return the number of nodes in the vessel
       @return unsigned
	 */
	unsigned GetNumberOfNodes();

	/**
       @return mVesselSegments.size()
	 */
	unsigned GetNumberOfSegments();
	
	/**
       @return mVesselSegmentLocations[index]
	 */
	boost::shared_ptr<CaVesselSegment<DIM> > GetSegment(unsigned index);

	/**
       @return mVesselSegments
	 */
	std::vector<boost::shared_ptr<CaVesselSegment<DIM> > > GetSegments();

	/**
       @return shared pointer to the first node of the first segment
	 */
	boost::shared_ptr<VascularNode<DIM> > GetStartNode();
	
    /**
     *  Return true if the vessel has data corresponding to the input key.
     *
     *  @return bool
     */
    bool HasDataKey(const std::string& rKey) const;

    /**
     *  Return whether the vessel is connected to another vessel.
     */
    bool IsConnectedTo(boost::shared_ptr<CaVessel<DIM> > pOtherVessel);
    
	/**
       Remove segments from the ends of a vessel
	 */
	void RemoveSegments(SegmentLocation::Value location);

    /**
     *  Add data of any type to the segment using the identifying key
     */
    template<typename T> void SetData(const std::string& rKey, T value);

	/**
	 *  Over-write the vessel's non-spatial DataContainer
	 *
	 *  This can be useful when copying data from an existing node.
	 */
	void SetDataContainer(const VasculatureData& rDataContainer);

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
       @return boost::shared_ptr<CaVessel<DIM> >
	 */
	boost::shared_ptr<CaVessel<DIM> > Shared();
};

#endif /* CAVESSEL_HPP_ */
