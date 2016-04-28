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

#ifndef VASCULARNODE_HPP_
#define VASCULARNODE_HPP_

#include <vector>
#include <string>
#include <map>
#include <boost/enable_shared_from_this.hpp>

#include "NodeFlowProperties.hpp"
#include "VasculatureData.hpp"
#include "UblasVectorInclude.hpp"
#include "ChastePoint.hpp"
#include "Cell.hpp"

/**
 *  Forward declaration to allow segments to manage adding and removing themselves from nodes.
 */
template<unsigned DIM>
class CaVesselSegment;

/**
 * This is a class for vascular nodes.
 *
 * Nodes are point locations along a vessel. They are useful for describing the end positions of
 * straight line vessel segments.
 */
template<unsigned DIM>
class VascularNode : public boost::enable_shared_from_this<VascularNode<DIM> >
{

    /**
     * Allow segments to manage adding and removing themselves from nodes.
     */
    friend class CaVesselSegment<DIM> ;

private:

    /**
     * Location of a node.
     */
    ChastePoint<DIM> mLocation;

    /**
     * Pointer to an associated Cell.
     */
    CellPtr mpCell;

    /**
     * Container for generic node data.
     */
    VasculatureData mDataContainer;

    /**
     * Id tag, can be useful for storing segment-node relationships in the VesselNetwork class.
     */
    unsigned mId;

    unsigned mTempId;

    /**
     * Label tag, can be useful for identifying input and output nodes.
     */
    std::string mLabel;

    /**
     * Collection of pointers to Vessel Segments connected to this node.
     */
    std::vector<boost::weak_ptr<CaVesselSegment<DIM> > > mVesselSegments;

    /**
     * Radius of the vessel at this node
     */
    double mRadius;

    /**
     * The flow properties for the node
     */
    boost::shared_ptr<NodeFlowProperties> mpNodeFlowProperties;

    /**
     * Is the vessel allowed to extend at this node
     */
    bool mIsMigrating;

public:

    /**
     * Constructor
     *
     * Create a node using a ChastePoint to specify the location
     *
     * @param rLocation the location of the node
     */
    VascularNode(const ChastePoint<DIM>& rLocation);

    /**
     * Constructor
     *
     * Create a node using xyz coordinates
     *
     * @param v1  the node's x-coordinate (defaults to 0)
     * @param v2  the node's y-coordinate (defaults to 0)
     * @param v3  the node's z-coordinate (defaults to 0)
     */
    VascularNode(double v1 = 0.0, double v2 = 0.0, double v3 = 0.0);

    /**
     * Constructor
     *
     * Create a node using ublas c_vector
     *
     * @param location the node's location (defaults to 0.0)
     */
    VascularNode(c_vector<double, DIM> location);

    /**
     * Copy constructor
     *
     * @param rExistingNode the node to copy from
     */
    VascularNode(const VascularNode<DIM>& rExistingNode);

    /**
     * Destructor
     */
    ~VascularNode();

    /**
     * Construct a new instance of the class and return a shared pointer to it.
     *
     * This method is included so that nodes can be created in a way that is consistent with
     * other vasculature features.
     *
     * @param rLocation the location of the node
     * @return a pointer to the newly created node
     */
    static boost::shared_ptr<VascularNode<DIM> > Create(const ChastePoint<DIM>& rLocation);

    /**
     * Construct a new instance of the class and return a shared pointer to it.
     *
     * This method is included so that nodes can be created in a way that is consistent with
     * other vasculature features.
     *
     * @param v1  the node's x-coordinate (defaults to 0)
     * @param v2  the node's y-coordinate (defaults to 0)
     * @param v3  the node's z-coordinate (defaults to 0)
     * @return a pointer to the newly created node
     */
    static boost::shared_ptr<VascularNode<DIM> > Create(double v1 = 0.0, double v2 = 0.0, double v3 = 0.0);

    /**
     * Construct a new instance of the class and return a shared pointer to it.
     *
     * This method is included so that nodes can be created in a way that is consistent with
     * other vasculature features.
     *
     * @param location the node's location (defaults to 0.0)
     * @return a pointer to the newly created node
     */
    static boost::shared_ptr<VascularNode<DIM> > Create(c_vector<double, DIM> location);

    /**
     * Construct a new instance of the class and return a shared pointer to it.
     *
     * This method is included so that nodes can be created in a way that is consistent with
     * other vasculature features.
     *
     * @param rExistingNode the node to copy from
     * @return a pointer to the newly created node
     */
    static boost::shared_ptr<VascularNode<DIM> > Create(const VascularNode<DIM>& rExistingNode);

    /**
     * Construct a new instance of the class and return a shared pointer to it.
     *
     * This method is included so that nodes can be created in a way that is consistent with
     * other vasculature features.
     *
     * @param pExistingNode the node to copy from
     * @return a pointer to the newly created node
     */
    static boost::shared_ptr<VascularNode<DIM> > Create(boost::shared_ptr<VascularNode<DIM> > pExistingNode);

    /**
     * Return a boost shared pointer to the associated Cell. Can be NULL if there is no cell assigned.
     *
     * @return a boost shared pointer to the associated Cell
     */
    CellPtr GetCell() const;

    /**
     * Return the node data for the input key. An attempt is made to cast to type T.
     *
     * @param rKey the key to be queried
     * @return the node data for the input key
     */
    template<typename T> T GetData(const std::string& rKey) const;

    /**
     * Return a const reference to the non-spatial data container.
     *
     * @return the data container for the node
     */
    const VasculatureData& rGetDataContainer() const;

    /**
     * Return a vector of data keys for the node. Input true if
     * the corresponding value should be castable to double.
     *
     * @param castableToDouble whether the returned keys should be castable to double
     * @return a vector of data keys for the node
     */
    std::vector<std::string> GetDataKeys(bool castableToDouble = false) const;

    /**
     * Return the distance between the input node and the node
     *
     * @param pNode a node to calculate the distance to
     * @return the distance to the node
     */
    double GetDistance(boost::shared_ptr<VascularNode<DIM> > pNode) const;

    /**
     * Return the distance between the input point and the node
     *
     * @param rPoint a ChastePoint at the location to calculate the distance to
     * @return the distance to the point
     */
    double GetDistance(const ChastePoint<DIM>& rPoint) const;

    /**
     * Return the distance between the input location and the node
     *
     * @param location the location to calculate the distance to
     * @return the distance to the location
     */
    double GetDistance(c_vector<double, DIM> location) const;

    /**
     * Return the flow properties of the node
     *
     * @return the flow properties of the node
     */
    boost::shared_ptr<NodeFlowProperties> GetFlowProperties() const;

    /**
     * Return the node Id
     *
     * @return the node id
     */
    unsigned GetId() const;

    unsigned GetTempId() const;

    /**
     * Return a const reference to the Label
     *
     * @return mLabel
     */
    const std::string& rGetLabel() const;

    /**
     * Return the location of the node or, if there is one, the associated Cell.
     *
     * @return a ChastePoint at the location of the node
     */
    ChastePoint<DIM> GetLocation() const;

    /**
     * Return the location of the node or, if there is one, the associated Cell.
     *
     * @return a ublas c_vector at the location of the node
     */
    c_vector<double, DIM> GetLocationVector() const;

    /**
     * Return the number of attached segments
     *
     * @return the number of segments attached to the node
     */
    unsigned GetNumberOfSegments() const;

    /**
     * Return the radius of the vessel at the node
     *
     * @return the radius of the vessel at the node
     */
    double GetRadius() const;

    /**
     * Return a map of nodal data for use by the vtk writer
     *
     * @return a map of nodal data for use by the vtk writer
     */
    std::map<std::string, double> GetVtkData() const;

    /**
     * Return a pointer to the specified vessel segment
     *
     * @param index the index of the vessel segment in the node's segment vector
     * @return a pointer to the segment for the input index
     */
    boost::shared_ptr<CaVesselSegment<DIM> > GetVesselSegment(unsigned index) const;

    /**
     * Return a vector of pointers to the attached vessel segments.
     *
     * @return a vector of pointers to the attached vessel segments
     */
    std::vector<boost::shared_ptr<CaVesselSegment<DIM> > > GetVesselSegments() const;

    /**
     * Return true if there is an associated Cell.
     *
     * @return whether there is a cell attached to the node.
     */
    bool HasCell() const;

    /**
     * Return true if the node has data corresponding to the input key.
     *
     * @param rKey the key to query the data collection with
     * @return whether the key exists in the data collection
     */
    bool HasDataKey(const std::string& rKey) const;

    /**
     * Return true if the input segment is attached to the node
     *
     * @param pSegment a pointer to the segment to query
     * @return whether the input segment is attached to the node
     */
    bool IsAttachedTo(const boost::shared_ptr<CaVesselSegment<DIM> > pSegment) const;

    /**
     * Return true if the node is coincident with the input location
     *
     * @param rPoint a ChastePoint at the query location
     * @return whether then node is coincident with the input location
     */
    bool IsCoincident(const ChastePoint<DIM>& rPoint) const;

    /**
     * Return true if the node is coincident with the input node
     *
     * @param pNode the node to test the location of
     * @return whether then node is coincident with the input node
     */
    bool IsCoincident(const boost::shared_ptr<VascularNode<DIM> > pNode) const;

    /**
     * Returns whether the node is actively migrating
     *
     * @return true if the node is actively migrating
     */
    bool IsMigrating() const;

    /**
     * Remove the assigned Cell
     */
    void RemoveCell();

    /**
     * Assign a Cell to the node. Overwrite any existing Cell.
     *
     * @param pCell the cell to assign to the node
     */
    void SetCell(CellPtr pCell);

    /**
     * Add data of any type to the node using the identifying key
     *
     * @param rKey the key for the data being assigned to the node
     */
    template<typename T> void SetData(const std::string& rKey, T value);

    /**
     * Over-write the node's non-spatial DataContainer
     *
     * This can be useful when copying data from an existing node.
     *
     * @param rDataContainer the data container to be copied from
     */
    void SetDataContainer(const VasculatureData& rDataContainer);

    /**
     * Assign the Id
     *
     * @param id the id for the node
     */

    void SetId(unsigned id);

    void SetTempId(unsigned id);

    /**
     * Assign the Label
     *
     * @param rLabel the label of the node
     */
    void SetLabel(const std::string& rLabel);

    /**
     * Set the flow properties of the node
     *
     * @param rFlowProperties the flow properties to be set
     */
    void SetFlowProperties(const NodeFlowProperties& rFlowProperties);

    /**
     * Set the vessel radius at this node
     *
     * @param radius the vessel radius
     */
    void SetRadius(double radius);

    /**
     * Set that the node is migrating
     * @param isMigrating whether the node is migrating
     */
    void SetIsMigrating(bool isMigrating);

    /**
     * Set the location of the node. This breaks any links with an assigned Cell, so if there is an
     * assigned Cell remove it.
     *
     * @param rLocation a ChastePoint at the location to be assigned
     */
    void SetLocation(const ChastePoint<DIM>& rLocation);

    /**
     * Set the location of the node. This breaks any links with an assigned Cell, so if there is an
     * assigned Cell remove it.
     *
     * @param location a ublas c_vector specifying the location
     */
    void SetLocation(c_vector<double, DIM> location);

private:

    /**
     * Add a vessel segment the node. Private because node-segment connectivity needs to be managed.
     *
     *  @param pVesselSegment the segment to be added
     */
    void AddSegment(boost::shared_ptr<CaVesselSegment<DIM> > pVesselSegment);

    /**
     * Remove a vessel segment from the node. Private because node-segment connectivity needs to be managed.
     *
     *  @param pVesselSegment the segment to be removed
     */
    void RemoveSegment(boost::shared_ptr<CaVesselSegment<DIM> > pVesselSegment);
};

#endif /* VASCULARNODE_HPP_ */
