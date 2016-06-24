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

#ifndef VesselNode_HPP_
#define VesselNode_HPP_

#include <vector>
#include <string>
#include <map>
#include <boost/enable_shared_from_this.hpp>
#include "AbstractVesselNetworkComponent.hpp"
#include "NodeFlowProperties.hpp"
#include "UblasVectorInclude.hpp"
#include "UnitCollection.hpp"
#include "VesselSegment.hpp"
#include "SmartPointers.hpp"

/**
 *  Forward declaration to allow segments to manage adding and removing themselves from nodes.
 */
template<unsigned DIM>
class VesselSegment;

/**
 * This is a class for vascular nodes, which are vessel network components.
 *
 * Nodes are point locations along a vessel. They are used for describing the end positions of
 * straight line vessel segments. Nodes are initialized in dimensionless units, however they can
 * be dimensionalized by setting a reference length scale.
 */
template<unsigned DIM>
class VesselNode : public boost::enable_shared_from_this<VesselNode<DIM> >, public AbstractVesselNetworkComponent<DIM>
{

    /**
     * Allow segments to manage adding and removing themselves from nodes.
     */
    friend class VesselSegment<DIM> ;

private:

    /**
     * Dimensionless location of a node in space.
     */
    ChastePoint<DIM> mLocation;

    /**
     * Collection of pointers to Vessel Segments connected to this node.
     */
    std::vector<boost::weak_ptr<VesselSegment<DIM> > > mSegments;

    /**
     * Is the vessel allowed to extend at this node
     */
    bool mIsMigrating;


    boost::shared_ptr<NodeFlowProperties<DIM> > mpFlowProperties;


    unsigned mPtrComparisonId;

    /**
     * Constructor. Kept private as the Create factory method should be used.
     *
     * Create a node using xyz coordinates
     *
     * @param v1  the node's x-coordinate (defaults to 0)
     * @param v2  the node's y-coordinate (defaults to 0)
     * @param v3  the node's z-coordinate (defaults to 0)
     */
    VesselNode(double v1 = 0.0, double v2 = 0.0, double v3 = 0.0);

    /**
     * Constructor. Kept private as the Create factory method should be used.
     *
     * Create a node using ublas c_vector
     *
     * @param location the node's location (defaults to 0.0)
     */
    VesselNode(c_vector<double, DIM> location);

    /**
     * Copy constructor. Kept private as the Create factory method should be used.
     *
     * @param rExistingNode the node to copy from
     */
    VesselNode(const VesselNode<DIM>& rExistingNode);

public:

    /**
     * Destructor
     */
    ~VesselNode();

    void SetComparisonId(unsigned id);

    unsigned GetComparisonId();

    /**
     * Construct a new instance of the class and return a shared pointer to it.
     *
     * This method is included so that nodes can be created in a way that is consistent with
     * other vasculature features
     *
     * @param v1  the node's x-coordinate (defaults to 0)
     * @param v2  the node's y-coordinate (defaults to 0)
     * @param v3  the node's z-coordinate (defaults to 0)
     * @return a pointer to the newly created node
     */
    static boost::shared_ptr<VesselNode<DIM> > Create(double v1 = 0.0, double v2 = 0.0, double v3 = 0.0);

    /**
     * Construct a new instance of the class and return a shared pointer to it.
     *
     * This method is included so that nodes can be created in a way that is consistent with
     * other vasculature features
     *
     * @param location the node's location (defaults to 0.0)
     * @return a pointer to the newly created node
     */
    static boost::shared_ptr<VesselNode<DIM> > Create(const c_vector<double, DIM>& location);

    /**
     * Construct a new instance of the class and return a shared pointer to it.
     *
     * This method is included so that nodes can be created in a way that is consistent with
     * other vasculature features
     *
     * @param rExistingNode the node to copy from
     * @return a pointer to the newly created node
     */
    static boost::shared_ptr<VesselNode<DIM> > Create(const VesselNode<DIM>& rExistingNode);

    /**
     * Construct a new instance of the class and return a shared pointer to it.
     *
     * This method is included so that nodes can be created in a way that is consistent with
     * other vasculature features.
     *
     * @param pExistingNode the node to copy from
     * @return a pointer to the newly created node
     */
    static boost::shared_ptr<VesselNode<DIM> > Create(boost::shared_ptr<VesselNode<DIM> > pExistingNode);

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
     * Return the flow properties of the component
     *
     * @return the flow properties of the component
     */
    boost::shared_ptr<NodeFlowProperties<DIM> > GetFlowProperties() const;

    /**
     * Return a reference to the scaled location of the node.
     *
     * @return a ublas c_vector at the location of the node
     */
    const c_vector<double, DIM>& rGetLocation() const;

    /**
     * Return a reference to the location of the node in SI units.
     *
     * @return a ublas c_vector at the location of the node
     */
    const c_vector<double, DIM>& rGetLocationSI() const;

    /**
     * Return the number of attached segments
     *
     * @return the number of segments attached to the node
     */
    unsigned GetNumberOfSegments() const;

    /**
     * Return a map of output data for writers
     * @return a map of component data for use by the vtk writer
     */
    std::map<std::string, double> GetOutputData();

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

    bool IsMigrating() const;

    /**
     * Set the flow properties of the node
     *
     * @param rFlowProperties the flow properties to be set
     */
    void SetFlowProperties(const NodeFlowProperties<DIM>& rFlowProperties);

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

#endif /* VesselNode_HPP_ */
