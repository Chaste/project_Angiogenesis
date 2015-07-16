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

#include <vector>
#include <string>
#include <boost/enable_shared_from_this.hpp>

#include "../../vessel/properties/SegmentFlowProperties.hpp"
#include "../../vessel/VasculatureData.hpp"
#include "UblasVectorInclude.hpp"
#include "ChastePoint.hpp"

/**
 *  Forward declaration to allow vessels to manage adding and removing themselves from segments.
 */
template<unsigned DIM>
class CaVessel;

template<unsigned DIM>
class VascularNode;

/**
 * This is a class for vessel segments.
 * .
 * Vessel segments are straight sub-units of vessels, defined by the positions of
 * their end nodes. Nodes cannot be created by the vessel segment class, they are
 * instead managed by the VascularNetwork class. Segments must always have two nodes.
 */
template<unsigned DIM>
class CaVesselSegment : public boost::enable_shared_from_this<CaVesselSegment<DIM> >
{
    /**
     * Allow vessels to manage adding and removing themselves from segments.
     */
    friend class CaVessel<DIM> ;

private:

    /**
     * Container for segment nodes
     */
    std::pair<boost::shared_ptr<VascularNode<DIM> >, boost::shared_ptr<VascularNode<DIM> > > mNodes;

    /**
     * Container for generic segment data.
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
     * Weak pointer to the vessel owning this segment
     */
    boost::weak_ptr<CaVessel<DIM> > mVessel;

    /**
     * Radius of the vessel at this segment
     */
    double mRadius;

    /**
     * The flow properties for the segment
     */
    boost::shared_ptr<SegmentFlowProperties> mpFlowProperties;

private:

    /**
     * Constructor - This is private as instances of this class must be created with a corresponding shared pointer. This is
     * implemented using the static Create method.
     *
     * @param pNode1 the first node in the segment
     * @param pNode2 the second node in the segment
     */
    CaVesselSegment(boost::shared_ptr<VascularNode<DIM> > pNode1, boost::shared_ptr<VascularNode<DIM> > pNode2);

public:

    /**
     * Copy Constructor - This should not be used directly as instances of this class must be created with a corresponding shared pointer. This is
     * implemented using the static Create method.
     *
     * @param rSegment the segment to be copied
     */
    CaVesselSegment(const CaVesselSegment<DIM>& rSegment);

    /**
     * Construct a new instance of the class and return a shared pointer to it. Also manage the association of segments to nodes by
     * passing self weak pointers to the nodes.
     *
     * @param pNode1 the first node in the segment
     * @param pNode2 the second node in the segment
     * @return a pointer to the newly created segment
     */
    static boost::shared_ptr<CaVesselSegment<DIM> > Create(boost::shared_ptr<VascularNode<DIM> > pNode1,
                                                           boost::shared_ptr<VascularNode<DIM> > pNode2);

    /**
     * Construct a new instance of the class and return a shared pointer to it. Also manage the association of segments to nodes by
     * passing self weak pointers to the nodes.
     *
     * @param pSegment the segment to be copied
     * @return a pointer to the newly created segment
     */
    static boost::shared_ptr<CaVesselSegment<DIM> > Create(boost::shared_ptr<CaVesselSegment<DIM> > pSegment);

    /*
     * Destructor
     */
    ~CaVesselSegment();

    /**
     * Copy a selection of member data and VasculatureData from the input segment. Convenient alternative to the copy constructor as nodes aren't
     * copied.
     * @param pTargetSegment the segment from which data is to be copied
     */
    void CopyDataFromExistingSegment(const boost::shared_ptr<CaVesselSegment<DIM> > pTargetSegment);

    void Remove();

    /**
     * Return the segment data for the input key. An attempt is made
     * to cast to type T.
     * @param rKey the key to be queried
     * @return the node data for the input key
     */
    template<typename T> T GetData(const std::string& rKey) const;

    /**
     * Return a const reference to the segment's non-spatial data container.
     *
     * @return the data container
     */
    const VasculatureData& rGetDataContainer() const;

    /**
     * Return a vector of data keys for the segment. Input true if
     * the corresponding value should be castable to double.
     *
     * @param castableToDouble whether the returned keys should be castable to double
     * @return a vector of data keys for the node
     */
    std::vector<std::string> GetDataKeys(bool castableToDouble = false) const;

    /**
     * Return the distance between the input point and the segment. If the projection of the
     * point is within the segment the distance is the perpendicular distance to the segment.
     * Otherwise it is the distance to the nearest node.
     *
     * @param rPoint the point the get the distance from
     * @return the distance to the segment
     */
    double GetDistance(const ChastePoint<DIM>& rPoint) const;

    /**
     * Return the distance between the input point and the segment. If the projection of the
     * point is within the segment the distance is the perpendicular distance to the segment.
     * Otherwise it is the distance to the nearest node.
     *
     * @param location the point the get the distance from
     * @return the distance to the segment
     */
    double GetDistance(c_vector<double, DIM> location) const;

    /**
     * Return the flow properties of the segment
     *
     * @return the flow properties of the segment
     */
    boost::shared_ptr<SegmentFlowProperties> GetFlowProperties() const;

    /**
     * Return the Id
     *
     * @return the segment id
     */
    unsigned GetId() const;

    /**
     * Return the Label
     *
     * @return the segment label
     */
    const std::string& rGetLabel() const;

    /**
     * Return the length
     *
     * @return the segment length
     */
    double GetLength() const;

    /**
     * Return the radius
     *
     * @return the segment radius
     */
    double GetRadius() const;

    /**
     * Return a point mid-way along the vessel segment
     *
     * @return a point midway along the segment
     */
    c_vector<double, DIM> GetMidPoint() const;

    /**
     * Return a pointer to the node specified by the index
     *
     * @return a pointer to the node specified by the index
     */
    boost::shared_ptr<VascularNode<DIM> > GetNode(unsigned index) const;

    /**
     * Return a pointer to the node on the other side of the segment
     *
     * @param pInputNode the node to get the opposite one to
     * @return a pointer to the node on the other side of the segment
     */
    boost::shared_ptr<VascularNode<DIM> > GetOppositeNode(boost::shared_ptr<VascularNode<DIM> > pInputNode) const;

    /**
     * Return the segment nodes as a pair
     *
     * @return the segment nodes as a pair
     */
    std::pair<boost::shared_ptr<VascularNode<DIM> >, boost::shared_ptr<VascularNode<DIM> > > GetNodes() const;

    /**
     * Return the projection of a point onto the segment. If the projection is outside the segment an
     * Exception is thrown.
     *
     * @param rPoint a ChastePoint at the location to be projected
     * @return the location of the projected point
     */
    c_vector<double, DIM> GetPointProjection(const ChastePoint<DIM>& rPoint) const;

    /**
     * Return the projection of a point onto the segment. If the projection is outside the segment an
     * Exception is thrown.
     *
     * @param location the location to be projected
     * @return the location of the projected point
     */
    c_vector<double, DIM> GetPointProjection(c_vector<double, DIM> location) const;

    /**
     * Return a unit vector pointing along the segment. The orientation along the segment is from node0 to node 1.
     *
     * @return a unit vector pointing along the segment
     */
    c_vector<double, DIM> GetUnitTangent() const;

    /**
     * Return a pointer to the vessel.
     *
     * @return the vessel attached to the segment
     */
    boost::shared_ptr<CaVessel<DIM> > GetVessel() const;

    /**
     * Return true if the segment has data corresponding to the input key.
     *
     * @param rKey the key for the data
     * @return whether the key exists in the data map
     */
    bool HasDataKey(const std::string& rKey) const;

    /**
     * Return whether the node is in the segment.
     *
     * @return whether the node is in the segment
     */
    bool HasNode(boost::shared_ptr<VascularNode<DIM> > pNode) const;

    /**
     * Return whether the segment is connected to another segment.
     *
     * @param pOtherSegment the segment to check connectivity with
     * @return whether the segment is connected to the input segment
     */
    bool IsConnectedTo(boost::shared_ptr<CaVesselSegment<DIM> > pOtherSegment) const;

    /**
     * Replace the node at the specified index with the passed in node.
     *
     * @param oldNodeIndex the index of the node to be replaced
     * @param pNewNode the node to be added to the segment
     */
    void ReplaceNode(unsigned oldNodeIndex, boost::shared_ptr<VascularNode<DIM> > pNewNode);

    /**
     *  Add data of any type to the segment using the identifying key
     *
     *  @param rKey the key associated with the data
     *  @param the value to be added to the data map
     */
    template<typename T> void SetData(const std::string& rKey, T value);

    /**
     *  Over-write the segment's non-spatial DataContainer
     *
     *  This can be useful when copying data from an existing segment.
     *
     *  @param rDataContainer the container to be inserted
     */
    void SetDataContainer(const VasculatureData& rDataContainer);

    /**
     * Set the flow properties of the segment
     *
     * @param rFlowProperties the flow properties to be set
     */
    void SetFlowProperties(const SegmentFlowProperties& rFlowProperties);

    /**
     * Assign the Id
     *
     * @param id the id to be assigned
     */
    void SetId(unsigned id);

    /**
     * Assign the Label
     *
     * @param rLabel the label to be assigned
     */
    void SetLabel(const std::string& rLabel);

    /**
     * Set the radius
     *
     * @radius the radius to be assigned
     */
    void SetRadius(double radius);

private:

    /**
     * Return a boost::shared_ptr to this object
     *
     * @return a shared pointer to the segment
     */
    boost::shared_ptr<CaVesselSegment<DIM> > Shared();

    /**
     * Add an adjoining Vessel to the segment.
     *
     * @param pVessel a vessel to be added to the segment
     */
    void AddVessel(boost::shared_ptr<CaVessel<DIM> > pVessel);

    /**
     * Remove an adjoining vessel from the segment.
     */
    void RemoveVessel();
};

#endif /* CAVESSELSEGMENT_HPP_ */
