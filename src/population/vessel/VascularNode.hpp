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
#include "UblasVectorInclude.hpp"
#include "UnitCollections.hpp"
#include "VesselSegment.hpp"
#include "SmartPointers.hpp"

/**
 *  Forward declaration to allow segments to manage adding and removing themselves from nodes.
 */
template<unsigned DIM>
class VesselSegment;

/**
 * This is a class for vascular nodes.
 *
 * Nodes are point locations along a vessel. They are useful for describing the end positions of
 * straight line vessel segments. Nodes are initialized in dimensionless units, however they can
 * be dimensionalized by setting a reference length scale.
 */
template<unsigned DIM>
class VascularNode : public boost::enable_shared_from_this<VascularNode<DIM> >
{

    /**
     * Allow segments to manage adding and removing themselves from nodes.
     */
    friend class VesselSegment<DIM> ;

private:

    /**
     * Location of a node. Metres implied.
     */
    ChastePoint<DIM> mLocation;

    /**
     * Container for generic node data.
     */
    std::map<std::string, double> mOutputData;

    /**
     * Id tag, useful for post-processing
     */
    unsigned mId;

    /**
     * Collection of pointers to Vessel Segments connected to this node.
     */
    std::vector<boost::weak_ptr<VesselSegment<DIM> > > mSegments;

    /**
     * Radius of the vessel at this node
     */
    units::quantity<unit::length> mRadius;

    /**
     * The flow properties for the node
     */
    boost::shared_ptr<NodeFlowProperties> mpFlowProperties;

    /**
     * Is the vessel allowed to extend at this node
     */
    bool mIsMigrating;

    /**
     * Reference length scale for dimensionalizing units
     */
    units::quantity<unit::length> mReferenceLength;

private:

    /**
     * Constructor. Kept private as the Create factory method should be used.
     *
     * Create a node using xyz coordinates, metres are assumed.
     *
     * @param v1  the node's x-coordinate (defaults to 0)
     * @param v2  the node's y-coordinate (defaults to 0)
     * @param v3  the node's z-coordinate (defaults to 0)
     */
    VascularNode(double v1 = 0.0, double v2 = 0.0, double v3 = 0.0);

    /**
     * Constructor. Kept private as the Create factory method should be used.
     *
     * Create a node using ublas c_vector, metres are assumed
     *
     * @param location the node's location (defaults to 0.0)
     */
    VascularNode(c_vector<double, DIM> location);

    /**
     * Copy constructor. Kept private as the Create factory method should be used.
     *
     * @param rExistingNode the node to copy from
     */
    VascularNode(const VascularNode<DIM>& rExistingNode);

public:

    /**
     * Destructor
     */
    ~VascularNode();

    /**
     * Construct a new instance of the class and return a shared pointer to it.
     *
     * This method is included so that nodes can be created in a way that is consistent with
     * other vasculature features, metres are assumed.
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
     * other vasculature features, metres are assumed.
     *
     * @param location the node's location (defaults to 0.0)
     * @return a pointer to the newly created node
     */
    static boost::shared_ptr<VascularNode<DIM> > Create(c_vector<double, DIM> location);

    /**
     * Construct a new instance of the class and return a shared pointer to it.
     *
     * This method is included so that nodes can be created in a way that is consistent with
     * other vasculature features
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
     * Return the non dimensional distance between the input location and the node
     *
     * @param rLocation the location to calculate the distance to
     * @return the distance to the location
     */
    double GetDistance(const c_vector<double, DIM>& rLocation) const;

    /**
     * Return the dimensional distance between the input location and the node
     *
     * @param rLocation the location to calculate the distance to
     * @return the distance to the location
     */
    units::quantity<unit::length> GetDimensionalDistance(const c_vector<double, DIM>& rLocation) const;

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

    /**
     * Return a refence to the dimensionless location of the node.
     *
     * @return a ublas c_vector at the location of the node
     */
    const c_vector<double, DIM>& rGetLocation() const;

    /**
     * Return a reference to the dimensional location of the node.
     *
     * @return a ublas c_vector at the location of the node
     */
    const c_vector<double, DIM>& rGetDimensionalLocation() const;

    /**
     * Return the number of attached segments
     *
     * @return the number of segments attached to the node
     */
    unsigned GetNumberOfSegments() const;

    /**
     * Return the output data for the given key.
     *
     * @param rKey the key to be queried
     * @return the node data for the input key
     */
    double GetOutputDataValue(const std::string& rKey) const;

    /**
     * Return a map of output data for writers
     * @return a map of nodal data for use by the vtk writer
     */
    std::map<std::string, double> GetOutputData() const;

    /**
     * Return the keys of the output data map
     * @param verbose include all flow data
     * @return a map of nodal data for use by the vtk writer
     */
    std::vector<std::string> GetOutputDataKeys() const;

    /**
     * Return the radius of the vessel at the node
     *
     * @return the radius of the vessel at the node
     */
    units::quantity<unit::length> GetDimensionalRadius() const;

    /**
     * Return the dimensionless radius of the vessel at the node
     *
     * @return the radius of the vessel at the node
     */
    double GetRadius() const;

    /**
     * Return the reference length for the node
     *
     * @return the reference length for the node
     */
    units::quantity<unit::length> GetReferenceLength() const;

    /**
     * Return the reference length for the node as a value and unit pair. This
     * incurs a cost relative to GetReferenceLength. It is used by the Python interface.
     *
     * @return a pair containing the reference length value and unit as a string.
     */
    std::pair<double, std::string> GetReferenceLengthValueAndUnit() const;

    /**
     * Return a pointer to the indexed vessel segment
     * @param index the segment index
     * @return a vector of pointers to the attached vessel segments
     */
    boost::shared_ptr<VesselSegment<DIM> > GetSegment(unsigned index) const;

    /**
     * Return a vector of pointers to the attached vessel segments.
     *
     * @return a vector of pointers to the attached vessel segments
     */
    std::vector<boost::shared_ptr<VesselSegment<DIM> > > GetSegments() const;

    /**
     * Return true if the input segment is attached to the node
     *
     * @param pSegment a pointer to the segment to query
     * @return whether the input segment is attached to the node
     */
    bool IsAttachedTo(const boost::shared_ptr<VesselSegment<DIM> > pSegment) const;

    /**
     * Return true if the node is coincident with the input location
     *
     * @param rLocation the query location
     * @return whether then node is coincident with the input location
     */
    bool IsCoincident(const c_vector<double, DIM>& rLocation) const;

    /**
     * Returns whether the node is actively migrating
     *
     * @return true if the node is actively migrating
     */
    bool IsMigrating() const;

    /**
     * Assign the Id
     *
     * @param id the id for the node
     */
    void SetId(unsigned id);

    /**
     * Set the flow properties of the node
     *
     * @param rFlowProperties the flow properties to be set
     */
    void SetFlowProperties(const NodeFlowProperties& rFlowProperties);

    /**
     * Set that the node is migrating
     * @param isMigrating whether the node is migrating
     */
    void SetIsMigrating(bool isMigrating);

    /**
     * Set the dimensionless location of the node.
     * @param rLocation a ublas c_vector specifying the location
     */
    void SetLocation(const c_vector<double, DIM>& rLocation);

    /**
     * Set the dimensionless location of the node.
     *
     * @param x the x location
     * @param y the y location
     * @param z the z location
     */
    void SetLocation(double x, double y, double z=0.0);

    /**
     * Add output data to the node using the identifying key
     *
     * @param rKey the key for the data being assigned to the node
     * @param value the value to be stored
     */
    void SetOutputData(const std::string& rKey, double value);

    /**
     * Set the dimensionless vessel radius at this node
     *
     * @param radius the vessel radius
     */
    void SetRadius(double radius);

    /**
     * Set the reference length for the node
     *
     * @param referenceLength the reference length
     */
    void SetReferenceLength(units::quantity<unit::length> referenceLength);

    /**
     * Set the reference length for the node. This method incurs a theoretical efficiency penalty relative to
     * the direct use of a Boost Unit. It is used by the Python interface.
     *
     * @param referenceLength the reference length
     * @param unit the length unit, currently "metres" and "microns" are supported.
     */
    void SetReferenceLength(double referenceLength, const std::string& unit);

    /**
     * The output data map can get out of date if radii or node indices change. It should be updated
     * before any write. This is not done in GetOutputData as that method needs to be const.
     */
    void UpdateOutputData();

private:

    /**
     * Add a vessel segment the node. Private because node-segment connectivity needs to be managed.
     *
     *  @param pVesselSegment the segment to be added
     */
    void AddSegment(boost::shared_ptr<VesselSegment<DIM> > pVesselSegment);

    /**
     * Remove a vessel segment from the node. Private because node-segment connectivity needs to be managed.
     *
     *  @param pVesselSegment the segment to be removed
     */
    void RemoveSegment(boost::shared_ptr<VesselSegment<DIM> > pVesselSegment);
};

#endif /* VASCULARNODE_HPP_ */
