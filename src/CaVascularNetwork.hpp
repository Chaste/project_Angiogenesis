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

#include <math.h>
#include <float.h>
#include <algorithm>
#include <utility>
#include <vector>
#include <iostream>
#define _BACKWARD_BACKWARD_WARNING_H 1 //Cut out the boost deprecated warning for now (gcc4.3)
#include <boost/config.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/multi_array.hpp>
#include <boost/enable_shared_from_this.hpp>
#include <boost/graph/visitors.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/pending/indirect_cmp.hpp>
#include "CaVessel.hpp"
#include "CaVascularNetworkNode.hpp"
#include "SmartPointers.hpp"

/*
    A vessel network is any set of connected or unconnected vessels within a model domain.
    The class contains a reference to a PottsMesh class which enables the vessels contained
    within a vessel network to interact spatially with other biological entities contained in the
    model domain.
 */

template<unsigned DIM>
class CaVascularNetwork : public boost::enable_shared_from_this<CaVascularNetwork<DIM> >
{

private:

    /**
       Container for Vessels in the VesselNetwork.
     */
    std::vector<boost::shared_ptr<CaVessel<DIM> > > mVesselArray;

    /**
       Container for VesselNetworkNodes in the VesselNetwork.
     */
    std::vector<boost::shared_ptr<CaVascularNetworkNode<DIM> > > mNodeArray;

    /**
       Level of haematocrit inside input vessels of a network.
     */
    double mArterialHaematocritLevel;

    /**
       Pressure at arterial input nodes of network. Pressure should be prescribed in units of
       Pascals.
     */
    double mArterialInputPressure;

    /**
       Pressure at venous output nodes of network. Pressure should be prescribed in units of
       Pascals.
     */
    double mVenousOutputPressure;

public:

    /*
     * Constructor
     *
     * Upon instantiation:
            The network contains no nodes or vessels.
            arterialHaematocritLevel is initialised as 0.45.
            arterialInputPressure is initialised as 3588 Pa.
            venousOutputPressure is initialised as 1993 Pa.
     */
    CaVascularNetwork();

    /*
     * Destructor
     */
    ~CaVascularNetwork();

    /**
       Returns a boost::shared_ptr to this VesselNetwork.

       @return boost::shared_ptr<CaVascularNetwork<DIM> >
     */
    boost::shared_ptr<CaVascularNetwork<DIM> > shared();

    /**
       This method is used to resolve the identity of vessels contained within a network.
       For example, a vessel could be fetched from the spatial mesh and provided as an argument to
       this function in order to determine the local identity of that vessel in the network.

       @param vessel boost::shared_ptr<CaVessel<DIM> >.
     */
    unsigned GetVesselID(boost::shared_ptr<CaVessel<DIM> > vessel);

    /**
       This method is used to resolve the identity of nodes contained within a network.
       For example, a node could be fetched from the spatial mesh and provided as an argument to
       this function in order to determine the local identity of that node in the network.

       @param node
       @return id of node (index of node in mNodeArray)
     */
    unsigned GetNodeID(boost::shared_ptr<CaVascularNetworkNode<DIM> > node);

    /**
       Returns the number of vessels contained within the network.
     */
    unsigned GetNumberOfVesselsInNetwork();

    /**
       Returns the number of nodes contained within the network.
     */
    unsigned GetNumberOfNodesInNetwork();

    /**
		Returns the number of vessels which occupy a particular location in the spatial mesh.
		At the moment, it is assumed that multiple vessels only occupy the same location on the
		spatial mesh if they are connected at that location. I.e. a node would also be present at
		that location.

		The spatial coordinate provided as an argument should be within the boundaries of the model
		domain (contained within the spatial mesh).

		@param location ChastePoint<DIM>.
     */
    unsigned GetNumberOfVesselsAtLocation(ChastePoint<DIM> coord);

    /**
		Returns a boost::shared_ptr to the vessel object with the prescribed id.

		@param vessel_id integer.
     */
    boost::shared_ptr<CaVessel<DIM> > GetVessel(int vessel_id);

    /**
		Returns a boost::shared_ptr to the vessel object with the prescribed location and at the
		prescribed position in the container at that location in the spatial mesh.

		@param location ChastePoint<DIM>.
		@param positionInContainer integer.
     */
    boost::shared_ptr<CaVessel<DIM> > GetVessel(ChastePoint<DIM> coord,
                                                 int positionInContainer);

    /**
		Returns a boost::shared_ptr to the node object with the prescribed id.

		@param node_id integer.
     */
    boost::shared_ptr<CaVascularNetworkNode<DIM> > GetNode(int node_id);

    /**
		Returns a boost::shared_ptr to the node object with the prescribed location in the spatial
		mesh. Currently, there should only ever be one node in any one location in the spatial mesh.

		To check whether there is a node at a location at all the method
		NumberOfNodesPresentAtLocation(location) may be used.
     */
    boost::shared_ptr<CaVascularNetworkNode<DIM> > GetNode(ChastePoint<DIM> location);

//    /**
//            Return mean length of vessels in network.
//     */
//    double GetMeanVesselLengthOfNeovasculature();
//
//    /**
//            Return mean radius of vessels in network.
//     */
//    double GetMeanVesselRadiusOfNeovasculature();
//
//    /**
//            Return mean tortuosity of vessels in network. Vessels whose tortuosity are infinite (i.e. self-loops)
//            are ignored in this calculation.
//     */
//    double GetMeanVesselTortuosityOfNeovasculature();
//
//    /**
//            Return number of vessels which have length within a specified range (lowerBoundLength < length <= upperBoundLength).
//     */
//    int GetNumberOfVesselsByLength(double lowerBoundLength, double upperBoundLength);
//
//    /**
//            Return number of vessels which have radius within a specified range (lowerBoundRadius < length <= upperBoundRadius).
//     */
//    int GetNumberOfVesselsByRadius(double lowerBoundRadius, double upperBoundRadius);
//
//    /**
//            Return number of vessels which have tortuosity within a specified range (lowerBoundTortuosity  < length <= upperBoundTortuosity).
//     */
//    int GetNumberOfVesselsByTortuosity(double lowerBoundTortuosity, double upperBoundTortuosity);

    /**
            Returns the haematocrit levels in all vessels connected to arterial input nodes.

            @return arterialHaematocritLevel.
     */
    double GetArterialHaematocritLevel();

    /**
            Returns the pressure at all input nodes in the network.

            @return arterialInputPressure.
     */
    double GetArterialInputPressure();

    /**
            Returns the pressure at all output nodes in the network.

            @return venousOutputPressure.
     */
    double GetVenousOutputPressure();

    /**
            Returns a array of vessels.

            @return VesselArray.
     */
    std::vector<boost::shared_ptr<CaVessel<DIM> > > GetVessels();

    /**
            Sets the member variable arterialHaematocritLevel to the prescribed value.
     */
    void SetArterialHaematocritLevel(double value);

    /**
            Sets the member variable arterialInputPressure to the prescribed value.
     */
    void SetArterialInputPressure(double value);

    /**
            Sets the member variable venousOutputPressure to the prescribed value.
     */
    void SetVenousOutputPressure(double value);

    /**
            Returns whether there is a node present at the prescribed location.
     */
    bool NodePresentAtLocation(ChastePoint<DIM> location);

    /**
            Returns how many nodes are present at a location. There should only ever really be at most
            one node present at any one location.  However, the VesselNetwork may enter a
            pseudo-acceptable state within some operation implementations whereby two nodes temporarily
            occupy the same location before being merged.
     */
    unsigned NumberOfNodesPresentAtLocation(ChastePoint<DIM> location);

    /**
            Adds a vessel to the VesselNetwork. A Vessel may only be added to the network if there is
            enough free space in the spatial mesh for the vessel to be added.

            The Vessel being added should not be adjoined to any other Vessel - joining to other vessels
            in the network occurs automatically within this method. Additionally, nodes defined within
            the Vessel being added are merged with existing nodes in the vessel network. Before adding
            the Vessel, therefore, the locations of the two nodes should be defined within that Vessel.

            This method also automatically adds the Vessel to the spatial mesh.
     */
    void AddVessel(boost::shared_ptr<CaVessel<DIM> > vessel);

    /**
            Checks whether the prescribed vessel is contained within the vessel network. This method is
            intended to be used mainly for debugging purposes.
     */
    bool VesselIsInNetwork(boost::shared_ptr<CaVessel<DIM> > vessel);

    /**
            Checks whether the prescribed node is contained within the vessel network. This method is
            intended to be used mainly for debugging purposes.
     */
    bool NodeIsInNetwork(boost::shared_ptr<CaVascularNetworkNode<DIM> > node);

    /**
            Return whether two nodes are connected.

            If the two nodes are the same this function returns true.
     */
    bool Connected(boost::shared_ptr<CaVascularNetworkNode<DIM> > node1, boost::shared_ptr<CaVascularNetworkNode<DIM> > node2);

    /**
            Return whether a node is connected to an input node.

            If the node is an input node this function returns true.
     */
    bool ConnectedToInputNode(boost::shared_ptr<CaVascularNetworkNode<DIM> > node);

    /**
            Return whether a nodes is connected to an output node.

            If the node is an output node this function returns true.
     */
    bool ConnectedToOutputNode(boost::shared_ptr<CaVascularNetworkNode<DIM> > node);

    /**
            Labels the node located at the prescribed location as an input node of the network.
            There must be a node located at the location provided.
            At present an input node must only be attached to one vessel.
     */
    void SetInputNode(ChastePoint<DIM> location);

    /**
            Labels the node located at the prescribed location as an output node of the network.
            There must be a node located at the location provided.
            At present an output node must only be attached to one vessel.
     */
    void SetOutputNode(ChastePoint<DIM> location);

    /**
            Save the VesselNetwork data to a file, identifiable by the prescribed string, in a format
            which may be imported to ParaView (file should be a .vtk file).
     */
    void SaveVasculatureDataToFile(string filename);
};

#endif /* CAVASCULARNETWORK_HPP_ */
