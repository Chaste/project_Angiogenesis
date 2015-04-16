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

	friend class CaVessel<DIM>;

private:

    /**
     * Container for segment end nodes
     */
    std::pair<boost::shared_ptr<VascularNode<DIM> >, boost::shared_ptr<VascularNode<DIM> > > mNodes;

    /**
     * Container for non-spatial segment data.
     */
	VasculatureData mDataContainer;

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

    /**
     *   Radius of the vessel at this segment
     */
    double mRadius;

    /**
     *   Haematocrit in the vessel at this segment
     */
    double mHaematocrit;

    /**
     *   Blood flow rate in the vessel at this segment
     */
    double mFlowRate;

    /**
     *   Impedance of this vessel segment
     */
    double mImpedance;


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
     *  Return the segment data for the input key. An attempt is made
     *  to cast to type T.
     *  @return T data
     */
    template<typename T> T GetData(const std::string& rKey);

	/**
	 *  Return a const reference to the segment's non-spatial data container.
	 *
	 *  @return mDataContainer
	 */
	const VasculatureData& rGetDataContainer() const;
	
    /**
     *  Return a vector of data keys for the segment. Input true if
     *  the corresponding value should be castable to double.
     *
     *  @return std::vector<std::string>
     */
    std::vector<std::string> GetDataKeys(bool castable_to_double = false) const;

    /**
     *  Return the distance between the input point and the segment. If the projection of the
     *  point is within the segment the distance is the perpendicular distance to the segment.
     *  Otherwise it is the distance to the nearest vascular node.
     *
     *  @return double
    */
    double GetDistance(const ChastePoint<DIM>& rPoint);

	/**
	 *  Return the Id
	 *
	 *  @return mId
	 */
	unsigned GetId() const;
	
	/**
	 *  Return the Label
	 *
	 *  @return mLabel
	 */
	const std::string& rGetLabel() const;

	/**
	 *  Return the length
	 */
    double GetLength() const;

	/**
	 *  Return the radius
	 */
    double GetRadius() const;

	/**
	 *  Return the haematocrit
	 */
    double GetHaematocrit() const;

	/**
	 *  Return the impedance
	 */
    double GetImpedance() const;

	/**
	 *  Return the flow rate
	 */
    double GetFlowRate() const;

	/**
	 *  Return a point mid-way along the vessel segment
	 */
    ChastePoint<DIM> GetMidPoint();
    
    /**
       Return a pointer to the node specified by the index
     */
    boost::shared_ptr<VascularNode<DIM> > GetNode(unsigned index);
    
    /**
       Return the segment nodes as a pair
     */
    std::pair<boost::shared_ptr<VascularNode<DIM> >, boost::shared_ptr<VascularNode<DIM> > > GetNodes();

    /**
     *  Return the projection of a point onto the segment. If the projection is outside the segment an
     *  Exception is thrown.
     */
    ChastePoint<DIM> GetPointProjection(const ChastePoint<DIM>& rPoint);

    /**
     *  Return a unit vector pointing along the segment. The orientation along the segment is from node0 to node 1.
     */
    ChastePoint<DIM> GetUnitTangent();

	/**
	 *  Return a boost::shared_ptr to the vessel.
	 *
	 * @return mAdjoiningVessels[i]
	 */
	boost::shared_ptr<CaVessel<DIM> > GetVessel();
	
    /**
     *  Return true if the segment has data corresponding to the input key.
     *
     *  @return bool
     */
    bool HasDataKey(const std::string& rKey) const;

    /**
     *  Return whether the node is in the segment.
     */
    bool HasNode(boost::shared_ptr<VascularNode<DIM> > pNode);

    /**
     *  Return whether the segment is connected to another segment.
     */
    bool IsConnectedTo(boost::shared_ptr<CaVesselSegment<DIM> > pOtherSegment);

    /**
       Replace the node at the specified index with the passed in node.
     */
    void ReplaceNode(unsigned oldNodeIndex, boost::shared_ptr<VascularNode<DIM> >  pNewNode);

    /**
     *  Add data of any type to the segment using the identifying key
     */
    template<typename T> void SetData(const std::string& rKey, T value);

	/**
	 *  Over-write the segment's non-spatial DataContainer
	 *
	 *  This can be useful when copying data from an existing segment.
	 */
	void SetDataContainer(const VasculatureData& rDataContainer);

	/**
	 *  Assign the Id
	 */
	void SetId(unsigned id);

	/**
	 *  Set the radius
	 */
    void SetRadius(double radius);

	/**
	 *  Set the haematocrit
	 */
    void SetHaematocrit(double haematocrit);

	/**
	 *  Set the flow rate
	 */
    void SetFlowRate(double flowRate);

	/**
	 *  Set the impedance
	 */
    void SetImpedance(double impedance);

	/**
	 *  Assign the Label
	 */
	void SetLabel(const std::string& rLabel);

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
