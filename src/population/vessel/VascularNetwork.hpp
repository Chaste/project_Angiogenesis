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

#ifndef VascularNetwork_HPP_
#define VascularNetwork_HPP_

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
#include "VascularNode.hpp"
#include "VasculatureData.hpp"
#include "UblasIncludes.hpp"

/**
 A vessel network is a collection of vessels
 */
template<unsigned DIM>
class VascularNetwork : public boost::enable_shared_from_this<VascularNetwork<DIM> >
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
    std::vector<boost::shared_ptr<VascularNode<DIM> > > mNodes;

    /**
     *  Is the data in mNodes up to date.
     */
    bool mNodesUpToDate;

    /**
     Container for vessel nodes in the VesselNetwork.
     */
    std::vector<boost::shared_ptr<VascularNode<DIM> > > mVesselNodes;

    /**
     *  Is the data in mVesselNodes up to date.
     */
    bool mVesselNodesUpToDate;

    /**
     * Container for non-spatial vessel network data.
     */
    VasculatureData mDataContainer;

public:

    /*
     * Constructor
     */
    VascularNetwork();

    /*
     * Construct a new instance of the class and return a shared pointer to it.
     */
    static boost::shared_ptr<VascularNetwork<DIM> > Create();

    /*
     * Destructor
     */
    ~VascularNetwork();

    /**
     Adds a vessel to the VesselNetwork.
     */
    void AddVessel(boost::shared_ptr<Vessel<DIM> > pVessel);

    /**
     Adds a collection of vessels to the VesselNetwork
     */
    void AddVessels(std::vector<boost::shared_ptr<Vessel<DIM> > > vessels);

    /**
     Make a copy of all vessels, but with new nodes and segments in each copy. Return the new vessels.
     */
    std::vector<boost::shared_ptr<Vessel<DIM> > > CopyVessels();

    /**
     Copy flow properties from the specified segment to all other segments
     */
    void CopySegmentFlowProperties(unsigned index=0);

    /**
     Make a copy of the selected vessels, but with new nodes and segments in each copy. Return the new vessels.
     */
    std::vector<boost::shared_ptr<Vessel<DIM> > > CopyVessels(std::vector<boost::shared_ptr<Vessel<DIM> > > vessels);

    /*
     * Divides a vessel into two at the specified location.
     */
    boost::shared_ptr<VascularNode<DIM> > DivideVessel(boost::shared_ptr<Vessel<DIM> > pVessel, ChastePoint<DIM> location);

    /*
     * Add a new node to the end of the vessel
     * @param pEndNode the node that the new segment will start on, should already be on the end of the vessel
     * @param pNewNode the new node to be added to the end of the vessel
     */
    void ExtendVessel(boost::shared_ptr<Vessel<DIM> > pVessel,
                      boost::shared_ptr<VascularNode<DIM> > pEndNode,
                      boost::shared_ptr<VascularNode<DIM> > pNewNode);

    /*
     * Forms a sprout at the specified locations.
     */
    boost::shared_ptr<Vessel<DIM> > FormSprout(ChastePoint<DIM> sproutBaseLocation, ChastePoint<DIM> sproutTipLocation);

    /**
     Get distance to nearest node
     */
    double GetDistanceToNearestNode(const ChastePoint<DIM>& rLocation);

    /**
     Get the node nearest to the specified location
     */
    boost::shared_ptr<VascularNode<DIM> > GetNearestNode(const ChastePoint<DIM>& rLocation);

    /**
     Get the node nearest to the specified location
     */
    boost::shared_ptr<VascularNode<DIM> > GetNearestNode(c_vector<double, DIM> location);

    /**
     Get the node nearest to the specified node
     */
    boost::shared_ptr<VascularNode<DIM> > GetNearestNode(boost::shared_ptr<VascularNode<DIM> > pInputNode);

    /**
     Get the segment nearest to the specified segment and the distance to it
     */
    std::pair<boost::shared_ptr<VesselSegment<DIM> >, double> GetNearestSegment(boost::shared_ptr<VesselSegment<DIM> > pSegment);

    /**
     Get the segment nearest to the specified node and the distance to it
     */
    std::pair<boost::shared_ptr<VesselSegment<DIM> >, double> GetNearestSegment(boost::shared_ptr<VascularNode<DIM> > pNode, bool sameVessel = true);

    /**
     Get the segment nearest to the specified location and the distance to it
     */
    std::pair<boost::shared_ptr<VesselSegment<DIM> >, double> GetNearestSegment(const ChastePoint<DIM>& rLocation);

    /**
     Get the segment nearest to the specified location and the distance to it
     */
    std::pair<boost::shared_ptr<VesselSegment<DIM> >, double> GetNearestSegment(c_vector<double, DIM> location);

    /**
     Get the intercapillary distance using a 2d measure
     */
    std::vector<double> GetInterCapillaryDistances();

    /**
     Get the segment nearest to the specified location
     */
    boost::shared_ptr<Vessel<DIM> > GetNearestVessel(const ChastePoint<DIM>& rLocation);

    /**
     Get the segment nearest to the specified location
     */
    boost::shared_ptr<Vessel<DIM> > GetNearestVessel(c_vector<double, DIM> location);

    /**
     Get the number of nodes near to a specified point
     */
    unsigned NumberOfNodesNearLocation(const ChastePoint<DIM>& rLocation, double radius = 0.0);

    /**
     Return the extents of the vessel network in the form ((xmin, xmax), (ymin, ymax), (zmin, zmax))
     */
    std::vector<std::pair<double, double> > GetExtents();

    /**
     Return the nodes in the network
     */
    std::vector<boost::shared_ptr<VascularNode<DIM> > > GetNodes();

    /**
     Return the node with the prescribed index
     */
    boost::shared_ptr<VascularNode<DIM> > GetNode(unsigned index);

    /**
     Return the node with the prescribed index
     */
    unsigned GetNodeIndex(boost::shared_ptr<VascularNode<DIM> > node);

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
    double GetTotalLength();

    /**
     Return the total volume of the network
     */
    double GetTotalVolume();

    /**
     Return the total surface area of the network
     */
    double GetTotalSurfaceArea();

    /**
     Return the average distance between segments
     */
    double GetAverageInterSegmentDistance();

    /**
     Return the average vessel length
     */
    double GetAverageVesselLength();

    /**
     Return a histogram of vessel length distributions
     */
    std::vector<unsigned> GetVesselLengthDistribution(double binSpacing = 10.0, unsigned numberOfBins = 10);

    /**
     Return the vessel with the specified index in the network
     */
    boost::shared_ptr<Vessel<DIM> > GetVessel(unsigned index);

    /**
     Return the only the nodes at the ends of vessels in the network
     */
    std::vector<boost::shared_ptr<VascularNode<DIM> > > GetVesselEndNodes();

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
    bool IsConnected(boost::shared_ptr<VascularNode<DIM> > pSourceNode,
                     boost::shared_ptr<VascularNode<DIM> > pQueryNode);

    /**
     Return whether a vector of nodes is connected to a vector of source nodes.
     */
    std::vector<bool> IsConnected(std::vector<boost::shared_ptr<VascularNode<DIM> > > sourceNodes,
                                  std::vector<boost::shared_ptr<VascularNode<DIM> > > queryNodes);

    /**
     * Return whether node is in network.
     */
    bool NodeIsInNetwork(boost::shared_ptr<VascularNode<DIM> > pSourceNode);

    /**
     Merge short vessels in the network
     */
    void MergeShortVessels(double cutoff = 10.0);

    /**
     Merge nodes with the same spatial location. Useful for
     tidying up networks read from file.
     */
    void MergeCoincidentNodes(double tolerance = 0.0);

    /**
     Merge nodes with the same spatial location. Useful for
     tidying up networks read from file.
     */
    void MergeCoincidentNodes(std::vector<boost::shared_ptr<Vessel<DIM> > > pVessels, double tolerance = 0.0);

    /**
     Merge nodes with the same spatial location. Useful for
     tidying up networks read from file.
     */
    void MergeCoincidentNodes(std::vector<boost::shared_ptr<VascularNode<DIM> > > nodes, double tolerance = 0.0);

    /*
     * Removes a vessel from the network
     * @param deleteVessel also remove the vessel from its child segments and nodes if true.
     */
    void RemoveVessel(boost::shared_ptr<Vessel<DIM> > pVessel, bool deleteVessel = false);

    /**
     Remove short vessels from the network
     */
    void RemoveShortVessels(double cutoff = 10.0, bool endsOnly = true);

    /**
     Apply the input data to all nodes in the network
     */
    void SetNodeData(VasculatureData data);

    /**
     Set the nodal radii to the same value
     */
    void SetNodeRadii(double radius);

    /**
     Set the properties of the segments in the network based on those of the prototype
     */
    void SetSegmentProperties(boost::shared_ptr<VesselSegment<DIM> > prototype);

    /**
     Apply the input data to all vessel segments in the network
     */
    void SetVesselData(VasculatureData data);

    /**
     Apply the input data to all vessels in the network
     */
    void SetSegmentData(VasculatureData data);

    /**
     Set the segment radii to the same value
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

#endif /* VascularNetwork_HPP_ */
