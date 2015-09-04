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

#ifndef CAVASCULARNETWORK_HPP_
#define CAVASCULARNETWORK_HPP_

#include <vector>
#include <set>
#include <map>
#ifdef CHASTE_VTK
#define _BACKWARD_BACKWARD_WARNING_H 1
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#endif // CHASTE_VTK
#include "CaVessel.hpp"
#include "CaVesselSegment.hpp"
#include "VascularNode.hpp"
#include "VasculatureData.hpp"
#include "UblasIncludes.hpp"

/**
 A vessel network is a collection of vessels
 */
template<unsigned DIM>
class CaVascularNetwork : public boost::enable_shared_from_this<CaVascularNetwork<DIM> >
{

private:

    /**
     Container for Vessels in the VesselNetwork.
     */
    std::vector<boost::shared_ptr<CaVessel<DIM> > > mVessels;

    /**
     Container for vessel segments in the VesselNetwork.
     */
    std::vector<boost::shared_ptr<CaVesselSegment<DIM> > > mSegments;

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
    CaVascularNetwork();

    /*
     * Construct a new instance of the class and return a shared pointer to it.
     */
    static boost::shared_ptr<CaVascularNetwork<DIM> > Create();

    /*
     * Destructor
     */
    ~CaVascularNetwork();

    /**
     Adds a vessel to the VesselNetwork.
     */
    void AddVessel(boost::shared_ptr<CaVessel<DIM> > pVessel);

    /**
     Adds a collection of vessels to the VesselNetwork
     */
    void AddVessels(std::vector<boost::shared_ptr<CaVessel<DIM> > > vessels);

    /**
     Make a copy of all vessels, but with new nodes and segments in each copy. Return the new vessels.
     */
    std::vector<boost::shared_ptr<CaVessel<DIM> > > CopyVessels();

    /**
     Make a copy of the selected vessels, but with new nodes and segments in each copy. Return the new vessels.
     */
    std::vector<boost::shared_ptr<CaVessel<DIM> > > CopyVessels(std::vector<boost::shared_ptr<CaVessel<DIM> > > vessels);

    /*
     * Removes a vessel from the network
     */
    void RemoveVessel(boost::shared_ptr<CaVessel<DIM> > pVessel, bool deleteVessel = false);

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
     Get the segment nearest to the specified segment and the distance to it
     */
    std::pair<boost::shared_ptr<CaVesselSegment<DIM> >, double> GetNearestSegment(boost::shared_ptr<CaVesselSegment<DIM> > pSegment);

    /**
     Get the segment nearest to the specified node and the distance to it
     */
    std::pair<boost::shared_ptr<CaVesselSegment<DIM> >, double> GetNearestSegment(boost::shared_ptr<VascularNode<DIM> > pNode);

    /**
     Get the segment nearest to the specified location and the distance to it
     */
    std::pair<boost::shared_ptr<CaVesselSegment<DIM> >, double> GetNearestSegment(const ChastePoint<DIM>& rLocation);

    /**
     Get the segment nearest to the specified location and the distance to it
     */
    std::pair<boost::shared_ptr<CaVesselSegment<DIM> >, double> GetNearestSegment(c_vector<double, DIM> location);

    std::vector<double> GetInterCapillaryDistances();

    /**
     Get the segment nearest to the specified location
     */
    boost::shared_ptr<CaVessel<DIM> > GetNearestVessel(const ChastePoint<DIM>& rLocation);

    /**
     Get the segment nearest to the specified location
     */
    boost::shared_ptr<CaVessel<DIM> > GetNearestVessel(c_vector<double, DIM> location);

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

    double GetTotalLength();

    double GetTotalVolume();

    double GetTotalSurfaceArea();

    double GetAverageInterSegmentDistance();

    double GetAverageVesselLength();

    std::vector<unsigned> GetVesselLengthDistribution(double binSpacing = 10.0, unsigned numberOfBins = 10);

    void RemoveShortVessels(double cutoff = 10.0, bool endsOnly = true);

    void MergeShortVessels(double cutoff = 10.0);

    /**
     Return the vessel with the specified index in the network
     */
    boost::shared_ptr<CaVessel<DIM> > GetVessel(unsigned index);

    /**
     Return the only the nodes at the ends of vessels in the network
     */
    std::vector<boost::shared_ptr<VascularNode<DIM> > > GetVesselEndNodes();

    /**
     Return the Index of the specified vessel
     */
    unsigned GetVesselIndex(boost::shared_ptr<CaVessel<DIM> > pVessel);

    /**
     Return the Index of the specified vessel segment
     */
    unsigned GetVesselSegmentIndex(boost::shared_ptr<CaVesselSegment<DIM> > pVesselSegment);

    /**
     Return the vessel segments in the network
     */
    std::vector<boost::shared_ptr<CaVesselSegment<DIM> > > GetVesselSegments();

    /**
     Return the vessels in the network
     */
    std::vector<boost::shared_ptr<CaVessel<DIM> > > GetVessels();

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
     Merge nodes with the same spatial location. Useful for
     tidying up networks read from file.
     */
    void MergeCoincidentNodes(double tolerance = 0.0);

    /**
     Merge nodes with the same spatial location. Useful for
     tidying up networks read from file.
     */
    void MergeCoincidentNodes(std::vector<boost::shared_ptr<CaVessel<DIM> > > pVessels, double tolerance = 0.0);

    /**
     Merge nodes with the same spatial location. Useful for
     tidying up networks read from file.
     */
    void MergeCoincidentNodes(std::vector<boost::shared_ptr<VascularNode<DIM> > > nodes, double tolerance = 0.0);

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
    void SetSegmentProperties(boost::shared_ptr<CaVesselSegment<DIM> > prototype);

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
    void Translate(const c_vector<double, DIM>& rTranslationVector, std::vector<boost::shared_ptr<CaVessel<DIM> > > vessels);

    /*
     * Divides a vessel into two at the specified location.
     */
    boost::shared_ptr<VascularNode<DIM> > DivideVessel(boost::shared_ptr<CaVessel<DIM> > pVessel, ChastePoint<DIM> location);

    /*
     * Forms a sprout at the specified locations.
     */
    boost::shared_ptr<CaVessel<DIM> > FormSprout(ChastePoint<DIM> sproutBaseLocation, ChastePoint<DIM> sproutTipLocation);

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

    /**
     Write the VesselNetwork data to a file.
     */
    void Write(const std::string& rFilename);

    /**
     * Outputs connectivity of vessels to file in graphviz format (.gv).
     */
    void WriteConnectivity(const std::string& rFilename);

};

#endif /* CAVASCULARNETWORK_HPP_ */
