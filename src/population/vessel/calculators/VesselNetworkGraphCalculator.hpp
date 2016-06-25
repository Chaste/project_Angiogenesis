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
#ifdef CHASTE_VTK
#define _BACKWARD_BACKWARD_WARNING_H 1
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#endif // CHASTE_VTK
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
     Container for Vessels in the VesselNetwork.
     */
    std::vector<boost::shared_ptr<Vessel<DIM> > > mVessels;

    /**
     Container for vessel segments in the VesselNetwork.
     */
    std::vector<boost::shared_ptr<VesselSegment<DIM> > > mSegments;

    /**
     *  Is the data in mSegments up to date.
     */
    bool mSegmentsUpToDate;

    /**
     Container for nodes in the VesselNetwork.
     */
    std::vector<boost::shared_ptr<VesselNode<DIM> > > mNodes;

    /**
     *  Is the data in mNodes up to date.
     */
    bool mNodesUpToDate;

    /**
     Container for vessel nodes in the VesselNetwork.
     */
    std::vector<boost::shared_ptr<VesselNode<DIM> > > mVesselNodes;

    /**
     *  Is the data in mVesselNodes up to date.
     */
    bool mVesselNodesUpToDate;

public:

    /*
     * Constructor
     */
    VesselNetwork();

    /*
     * Construct a new instance of the class and return a shared pointer to it.
     */
    static boost::shared_ptr<VesselNetwork<DIM> > Create();

    /*
     * Destructor
     */
    ~VesselNetwork();

    /**
     Adds a vessel to the VesselNetwork.
     */
    void AddVessel(boost::shared_ptr<Vessel<DIM> > pVessel);

    /**
     Adds a collection of vessels to the VesselNetwork
     */
    void AddVessels(std::vector<boost::shared_ptr<Vessel<DIM> > > vessels);

    /**
     Copy flow properties from the specified segment to all other segments
     */
    void CopySegmentFlowProperties(unsigned index=0);

    /**
     Make a copy of all vessels, but with new nodes and segments in each copy. Return the new vessels.
     */
    std::vector<boost::shared_ptr<Vessel<DIM> > > CopyVessels();

    /**
     Make a copy of the selected vessels, but with new nodes and segments in each copy. Return the new vessels.
     */
    std::vector<boost::shared_ptr<Vessel<DIM> > > CopyVessels(std::vector<boost::shared_ptr<Vessel<DIM> > > vessels);

    /*
     * Divides a vessel into two at the specified location.
     */
    boost::shared_ptr<VesselNode<DIM> > DivideVessel(boost::shared_ptr<Vessel<DIM> > pVessel, const c_vector<double, DIM>& location);

    /*
     * Add a new node to the end of the vessel
     * @param pEndNode the node that the new segment will start on, should already be on the end of the vessel
     * @param pNewNode the new node to be added to the end of the vessel
     */
    void ExtendVessel(boost::shared_ptr<Vessel<DIM> > pVessel, boost::shared_ptr<VesselNode<DIM> > pEndNode,
                      boost::shared_ptr<VesselNode<DIM> > pNewNode);

    /*
     * Forms a sprout at the specified locations.
     */
    boost::shared_ptr<Vessel<DIM> > FormSprout(const c_vector<double, DIM>& sproutBaseLocation,
                                               const c_vector<double, DIM>& sproutTipLocation);

    /**
     Get distance to nearest node
     */
    double GetDistanceToNearestNode(const c_vector<double, DIM>& rLocation);

    /**
     Get the node nearest to the specified location
     */
    boost::shared_ptr<VesselNode<DIM> > GetNearestNode(const c_vector<double, DIM>& rLocation);

    /**
     Get the node nearest to the specified node
     */
    boost::shared_ptr<VesselNode<DIM> > GetNearestNode(boost::shared_ptr<VesselNode<DIM> > pInputNode);

    /**
     Get the segment nearest to the specified segment and the distance to it
     */
    std::pair<boost::shared_ptr<VesselSegment<DIM> >, double> GetNearestSegment(boost::shared_ptr<VesselSegment<DIM> > pSegment);

    /**
     Get the segment nearest to the specified node and the distance to it
     */
    std::pair<boost::shared_ptr<VesselSegment<DIM> >, double> GetNearestSegment(boost::shared_ptr<VesselNode<DIM> > pNode, bool sameVessel = true);

    /**
     Get the segment nearest to the specified location and the distance to it
     */
    std::pair<boost::shared_ptr<VesselSegment<DIM> >, double> GetNearestSegment(const c_vector<double, DIM>& location);

    /**
     Get the intercapillary distance using a 2d measure
     */
    std::vector<units::quantity<unit::length> > GetInterCapillaryDistances();

    /**
     Get the segment nearest to the specified location
     */
    boost::shared_ptr<Vessel<DIM> > GetNearestVessel(const c_vector<double, DIM>& location);

    /**
     Get index of nearest node
     */
    unsigned GetNodeIndex(boost::shared_ptr<VesselNode<DIM> > node);

    /**
     Get the number of nodes near to a specified point
     */
    unsigned NumberOfNodesNearLocation(const c_vector<double, DIM>&  rLocation, double radius = 0.0 * unit::metres);

    /**
     Return the extents of the vessel network in the form ((xmin, xmax), (ymin, ymax), (zmin, zmax))
     */
    std::vector<std::pair<double, double> > GetExtents();

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
     Return the total length of the network
     */
    units::quantity<unit::length> GetTotalLength();

    /**
     Return the total volume of the network
     */
    units::quantity<unit::volume> GetTotalVolume();

    /**
     Return the total surface area of the network
     */
    units::quantity<unit::area> GetTotalSurfaceArea();

    /**
     Return the average distance between segments
     */
    units::quantity<unit::length> GetAverageInterSegmentDistance();

    /**
     Return the average vessel length
     */
    units::quantity<unit::length> GetAverageVesselLength();

    /**
     Return a histogram of vessel length distributions
     */
    std::vector<unsigned> GetVesselLengthDistribution(double binSpacing = 10.0, unsigned numberOfBins = 10);

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
     Return the vessels in the network
     */
    std::vector<boost::shared_ptr<Vessel<DIM> > > GetVessels();

    /**
     Return the indices of each node attached to a node
     */
    std::vector<std::vector<unsigned> > GetNodeNodeConnectivity();

    /**
     Return the indices of each vessel attached to a node
     */
    std::vector<std::vector<unsigned> > GetNodeVesselConnectivity();

    /**
     Return the the vessel network in vtk form
     */
    vtkSmartPointer<vtkPolyData> GetVtk();

    /**
     Return whether a node is connected to a source node.
     */
    bool IsConnected(boost::shared_ptr<VesselNode<DIM> > pSourceNode, boost::shared_ptr<VesselNode<DIM> > pQueryNode);

    /**
     Return whether a vector of nodes is connected to a vector of source nodes.
     */
    std::vector<bool> IsConnected(std::vector<boost::shared_ptr<VesselNode<DIM> > > sourceNodes,
                                  std::vector<boost::shared_ptr<VesselNode<DIM> > > queryNodes);

    /**
     * Return whether node is in network.
     */
    bool NodeIsInNetwork(boost::shared_ptr<VesselNode<DIM> > pSourceNode);

    /**
     * Merge short vessels in the network
     */
    void MergeShortVessels(double cutoff = 10.0e-6);

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

    /*
     * Removes a vessel from the network
     * @param deleteVessel also remove the vessel from its child segments and nodes if true.
     */
    void RemoveVessel(boost::shared_ptr<Vessel<DIM> > pVessel, bool deleteVessel = false);

    /**
     * Remove short vessels from the network
     *
     */
    void RemoveShortVessels(double cutoff = 10.0, bool endsOnly = true);

    /**
     * Set the nodal radii to the same value
     */
    void SetNodeRadii(double radius);

    /**
     * Set the properties of the segments in the network based on those of the prototype
     */
    void SetSegmentProperties(boost::shared_ptr<VesselSegment<DIM> > prototype);

    /**
     * Set the segment radii to the same value
     */
    void SetSegmentRadii(double radius);

    /*
     * Translate the network along the provided vector
     */
    void Translate(const c_vector<double, DIM>& rTranslationVector);

    /*
     * Translate specific vessels along the provided vector
     */
    void Translate(const c_vector<double, DIM>& rTranslationVector, std::vector<boost::shared_ptr<Vessel<DIM> > > vessels);

    /*
     * Update the network node collection
     */
    void UpdateNodes();

    /*
     * Update the network segment collection
     */
    void UpdateSegments();

    /*
     * Update the network vessel collection
     */
    void UpdateVesselNodes();

    /*
     * Update the vessel id tags
     */
    void UpdateVesselIds();

    /*
     * Update all dynamic storage in the vessel network, optionally merge coinciden nodes
     */
    void UpdateAll(bool merge=false);

    /**
     Write the VesselNetwork data to a file.
     */
    void Write(const std::string& rFilename);

    /**
     * Outputs connectivity of vessels to file in graphviz format (.gv).
     */
    void WriteConnectivity(const std::string& rFilename);

    /**
     * Returns whether a vessel crosses a line segment.
     */
    bool VesselCrossesLineSegment(c_vector<double, DIM> coordinate_1, c_vector<double, DIM> coordinate_2, double radius = 1e-6);

};

#endif /* VesselNetwork_HPP_ */
