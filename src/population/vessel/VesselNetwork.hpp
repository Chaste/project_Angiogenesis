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

#ifndef VesselNetwork_HPP_
#define VesselNetwork_HPP_

#include <vector>
#include <set>
#include <map>
#include "Vessel.hpp"
#include "VesselSegment.hpp"
#include "VesselNode.hpp"
#include "UblasIncludes.hpp"
#include "UnitCollection.hpp"
#include "AbstractVesselNetworkComponent.hpp"

/**
 * A vessel network is a collection of vessels.
 */
template<unsigned DIM>
class VesselNetwork : public boost::enable_shared_from_this<VesselNetwork<DIM> >, public AbstractVesselNetworkComponent<DIM>
{

private:

    /**
     * Container for Vessels in the VesselNetwork.
     */
    std::vector<boost::shared_ptr<Vessel<DIM> > > mVessels;

    /**
     * Container for vessel segments in the VesselNetwork.
     */
    std::vector<boost::shared_ptr<VesselSegment<DIM> > > mSegments;

    /**
     * Is the data in mSegments up to date.
     */
    bool mSegmentsUpToDate;

    /**
     * Container for nodes in the VesselNetwork.
     */
    std::vector<boost::shared_ptr<VesselNode<DIM> > > mNodes;

    /**
     * Is the data in mNodes up to date.
     */
    bool mNodesUpToDate;

    /**
     * Container for vessel nodes in the VesselNetwork.
     */
    std::vector<boost::shared_ptr<VesselNode<DIM> > > mVesselNodes;

    /**
     * Is the data in mVesselNodes up to date.
     */
    bool mVesselNodesUpToDate;

public:

    /**
     * Constructor.
     */
    VesselNetwork();

    /**
     * Construct a new instance of the class and return a shared pointer to it.
     */
    static boost::shared_ptr<VesselNetwork<DIM> > Create();

    /**
     * Destructor
     */
    ~VesselNetwork();

    /**
     * Adds a vessel to the VesselNetwork.
     */
    void AddVessel(boost::shared_ptr<Vessel<DIM> > pVessel);

    /**
     * Adds a collection of vessels to the VesselNetwork
     */
    void AddVessels(std::vector<boost::shared_ptr<Vessel<DIM> > > vessels);

    /**
     * Copy flow properties from the specified segment to all other segments
     */
    void CopySegmentFlowProperties(unsigned index=0);

    /**
     * Make a copy of all vessels, but with new nodes and segments in each copy. Return the new vessels.
     */
    std::vector<boost::shared_ptr<Vessel<DIM> > > CopyVessels();

    /**
     * Make a copy of the selected vessels, but with new nodes and segments in each copy. Return the new vessels.
     */
    std::vector<boost::shared_ptr<Vessel<DIM> > > CopyVessels(std::vector<boost::shared_ptr<Vessel<DIM> > > vessels);

    /**
     * Divides a vessel into two at the specified location.
     */
    boost::shared_ptr<VesselNode<DIM> > DivideVessel(boost::shared_ptr<Vessel<DIM> > pVessel,
                                                     const DimensionalChastePoint<DIM>& location);

    /**
     * Add a new node to the end of the vessel
     * @param pEndNode the node that the new segment will start on, should already be on the end of the vessel
     * @param pNewNode the new node to be added to the end of the vessel
     */
    void ExtendVessel(boost::shared_ptr<Vessel<DIM> > pVessel, boost::shared_ptr<VesselNode<DIM> > pEndNode,
                      boost::shared_ptr<VesselNode<DIM> > pNewNode);

    /**
     * Forms a sprout at the specified locations.
     */
    boost::shared_ptr<Vessel<DIM> > FormSprout(const DimensionalChastePoint<DIM>& sproutBaseLocation,
                                               const DimensionalChastePoint<DIM>& sproutTipLocation);

    /**
     * Get distance to nearest node
     */
    units::quantity<unit::length> GetDistanceToNearestNode(const DimensionalChastePoint<DIM>& rLocation);

    /**
     * Get the node nearest to the specified location
     */
    boost::shared_ptr<VesselNode<DIM> > GetNearestNode(const DimensionalChastePoint<DIM>& rLocation);

    /**
     Get the node nearest to the specified node
     */
    boost::shared_ptr<VesselNode<DIM> > GetNearestNode(boost::shared_ptr<VesselNode<DIM> > pInputNode);

    /**
     Get the segment nearest to the specified segment and the distance to it
     */
    std::pair<boost::shared_ptr<VesselSegment<DIM> >, units::quantity<unit::length> > GetNearestSegment(boost::shared_ptr<VesselSegment<DIM> > pSegment);

    /**
     Get the segment nearest to the specified node and the distance to it
     */
    std::pair<boost::shared_ptr<VesselSegment<DIM> >, units::quantity<unit::length> > GetNearestSegment(boost::shared_ptr<VesselNode<DIM> > pNode, bool sameVessel = true);

    /**
     Get the segment nearest to the specified location and the distance to it
     */
    std::pair<boost::shared_ptr<VesselSegment<DIM> >, units::quantity<unit::length> > GetNearestSegment(const DimensionalChastePoint<DIM>& location);

    /**
     Get the segment nearest to the specified location
     */
    boost::shared_ptr<Vessel<DIM> > GetNearestVessel(const DimensionalChastePoint<DIM>& location);

    /**
     Get index of nearest node
     */
    unsigned GetNodeIndex(boost::shared_ptr<VesselNode<DIM> > node);

    /**
     Get the number of nodes near to a specified point
     */
    unsigned NumberOfNodesNearLocation(const DimensionalChastePoint<DIM>&  rLocation, double tolerance = 0.0);

    /**
     Return the extents of the vessel network in the form ((xmin, xmax), (ymin, ymax), (zmin, zmax))
     */
    std::pair<DimensionalChastePoint<DIM>, DimensionalChastePoint<DIM> > GetExtents(bool useRadii = false);

    /**
     * Return the indexed node in the network
     * This is dangerous as the node member array can be updated and this index will be out of date, use GetNodes instead.
     */
    boost::shared_ptr<VesselNode<DIM> > GetNode(unsigned index);

    /**
     Return the nodes in the network
     */
    std::vector<boost::shared_ptr<VesselNode<DIM> > > GetNodes();

    /**
     Return the number of nodes in the network.
     */
    unsigned GetNumberOfNodes();

    /**
     Return the number of vessel nodes in the network.
     */
    unsigned GetNumberOfVesselNodes();

    /**
     Return the number of vessels in the network.
     */
    unsigned GetNumberOfVessels();

    /**
     Return the number of branches on the most highly connected node
     */
    unsigned GetMaxBranchesOnNode();

    /**
     Return the only the nodes at the ends of vessels in the network
     */
    std::vector<boost::shared_ptr<VesselNode<DIM> > > GetVesselEndNodes();

    /**
     Return the Index of the specified vessel
     */
    unsigned GetVesselIndex(boost::shared_ptr<Vessel<DIM> > pVessel);

    /**
     Return the Index of the specified vessel segment
     */
    unsigned GetVesselSegmentIndex(boost::shared_ptr<VesselSegment<DIM> > pVesselSegment);

    /**
     Return the vessel segments in the network
     */
    std::vector<boost::shared_ptr<VesselSegment<DIM> > > GetVesselSegments();

    /**
     Return the indexed vessel
     */
    boost::shared_ptr<Vessel<DIM> > GetVessel(unsigned index);

    /**
     Return the vessels in the network
     */
    std::vector<boost::shared_ptr<Vessel<DIM> > > GetVessels();

    /**
     * Return whether node is in network.
     */
    bool NodeIsInNetwork(boost::shared_ptr<VesselNode<DIM> > pSourceNode);

    /**
     * Merge short vessels in the network
     */
    void MergeShortVessels(units::quantity<unit::length> cutoff = 10.0 * 1.e-6 * unit::metres);

    /**
     * Merge nodes with the same spatial location. Useful for
     * tidying up networks read from file.
     */
    void MergeCoincidentNodes(double tolerance = 0.0);

    /**
     * Merge nodes with the same spatial location. Useful for
     * tidying up networks read from file.
     */
    void MergeCoincidentNodes(std::vector<boost::shared_ptr<Vessel<DIM> > > pVessels, double tolerance = 0.0);

    /**
     * Merge nodes with the same spatial location. Useful for
     * tidying up networks read from file.
     */
    void MergeCoincidentNodes(std::vector<boost::shared_ptr<VesselNode<DIM> > > nodes, double tolerance = 0.0);

    /**
     * Removes a vessel from the network
     * @param deleteVessel also remove the vessel from its child segments and nodes if true.
     */
    void RemoveVessel(boost::shared_ptr<Vessel<DIM> > pVessel, bool deleteVessel = false);

    /**
     * Remove short vessels from the network
     *
     */
    void RemoveShortVessels(units::quantity<unit::length> cutoff = 10.0* 1.e-6 * unit::metres, bool endsOnly = true);

    /**
     * Set the nodal radii to the same value
     */
    void SetNodeRadii(units::quantity<unit::length> radius);

    /**
     * Get the node radius by averaging its segments
     */
    void SetNodeRadiiFromSegments();

    /**
     * Set the properties of the segments in the network based on those of the prototype
     */
    void SetSegmentProperties(boost::shared_ptr<VesselSegment<DIM> > prototype);

    /**
     * Set the segment radii to the same value
     */
    void SetSegmentRadii(units::quantity<unit::length> radius);

    /**
     * Set the segment viscosity to the same value
     */
    void SetSegmentViscosity(units::quantity<unit::dynamic_viscosity> viscosity);

    /**
     * Translate the network along the provided vector
     */
    void Translate(DimensionalChastePoint<DIM> rTranslationVector);

    /**
     * Translate specific vessels along the provided vector
     */
    void Translate(DimensionalChastePoint<DIM> rTranslationVector, std::vector<boost::shared_ptr<Vessel<DIM> > > vessels);

    /**
     * Update the network node collection
     */
    void UpdateNodes();

    /**
     * Update the network segment collection
     */
    void UpdateSegments();

    /**
     * Update the network vessel collection
     */
    void UpdateVesselNodes();

    /**
     * Update the vessel id tags
     */
    void UpdateVesselIds();

    /**
     * Update all dynamic storage in the vessel network, optionally merge coinciden nodes
     */
    void UpdateAll(bool merge=false);

    /**
     * Returns whether a vessel crosses a line segment.
     */
    bool VesselCrossesLineSegment(DimensionalChastePoint<DIM> coordinate_1, DimensionalChastePoint<DIM> coordinate_2, double tolerance = 1e-6);

    /**
     * Write the network to file
     */
    void Write(const std::string& rFileName);

};

#endif /* VesselNetwork_HPP_ */
