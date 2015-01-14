/*
 * CaVascularNetwork.hpp
 *
 *  Created on: 13 Jan 2015
 *      Author: connor
 */

#ifndef CAVASCULARNETWORK_HPP_
#define CAVASCULARNETWORK_HPP_

#include "CaVessel.hpp"
#include "CaVascularNetworkNode.hpp"
#include "SpatialCoordinate.hpp"
#include "CaBasedCellPopulation.hpp"

#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>

#include <boost/shared_ptr.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/multi_array.hpp>
#include <boost/enable_shared_from_this.hpp>

//! CaVascularNetwork class

/*!
    A vessel network is any set of connected or unconnected vessels within a model domain.
    The class contains a reference to a PottsMesh class which enables the vessels contained
    within a vessel network to interact spatially with other biological entities contained in the
    model domain.
 */

template<unsigned SPATIAL_DIM>
class CaVascularNetwork : public boost::enable_shared_from_this<CaVascularNetwork<SPATIAL_DIM> >
{

protected:

    /**
       Container for Vessels contained within the VesselNetwork.
     */
    std::vector<boost::shared_ptr<CaVessel<SPATIAL_DIM> > > mVesselArray;

    /**
       Container for VesselNetworkNodes contained within the VesselNetwork.
     */
    std::vector<boost::shared_ptr<CaVascularNetworkNode<SPATIAL_DIM> > > mNodeArray;

    /**
       Level of haematocrit inside input vessels of a network.
     */
    double mArterialHaematocritLevel;

    /**
       Pressure at arterial input nodes of network. Pressure should be prescribed in units of
       Pascals.

       Nodes are prescribed as input nodes using GetNode(...)->SetIsInputNode(true) or
       SetInputNode(location)
     */
    double mArterialInputPressure;

    /**
       Pressure at venous output nodes of network. Pressure should be prescribed in units of
       Pascals.

       Nodes are prescribed as output nodes using GetNode(...)->SetIsOutputNode(true) or
       SetOutputNode(location)
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
    virtual ~CaVascularNetwork();

    /**
       Returns a boost::shared_ptr to this VesselNetwork.

       @return boost::shared_ptr<CaVascularNetwork<SPATIAL_DIM> >
     */
    boost::shared_ptr<CaVascularNetwork<SPATIAL_DIM> > shared();

    /**
       This method is used to resolve the identity of vessels contained within a network.
       For example, a vessel could be fetched from the spatial mesh and provided as an argument to
       this function in order to determine the local identity of that vessel in the network.

       @param vessel boost::shared_ptr<CaVessel<SPATIAL_DIM> >.
     */
    int GetVesselID(boost::shared_ptr<CaVessel<SPATIAL_DIM> > vessel);

    /**
       This method is used to resolve the identity of nodes contained within a network.
       For example, a node could be fetched from the spatial mesh and provided as an argument to
       this function in order to determine the local identity of that node in the network.

       @param node
       @return id of node (index of node in mNodeArray)
     */
    int GetNodeID(boost::shared_ptr<CaVascularNetworkNode<SPATIAL_DIM> > node);

    /**
       Returns the number of vessels contained within the network.
     */
    int GetNumberOfVesselsInNetwork();

    /**
            Returns the number of nodes contained within the network.
     */
    int GetNumberOfNodesInNetwork();

    /**
            Returns the number of vessels which occupy a particular location in the spatial mesh.
            At the moment, it is assumed that multiple vessels only occupy the same location on the
            spatial mesh if they are connected at that location. I.e. a node would also be present at
            that location.

            The spatial coordinate provided as an argument should be within the boundaries of the model
            domain (contained within the spatial mesh).

            @param location ChastePoint<SPATIAL_DIM>.
     */
    int GetNumberOfVesselsAtLocation(ChastePoint<SPATIAL_DIM> coord);

    /**
            Returns a boost::shared_ptr to the vessel object with the prescribed id.

            @param vessel_id integer.
     */
    boost::shared_ptr<CaVessel<SPATIAL_DIM> > GetVessel(int vessel_id);

    /**
            Returns a boost::shared_ptr to the vessel object with the prescribed location and at the
            prescribed position in the container at that location in the spatial mesh.

            @param location ChastePoint<SPATIAL_DIM>.
            @param positionInContainer integer.
     */
    boost::shared_ptr<CaVessel<SPATIAL_DIM> > GetVessel(ChastePoint<SPATIAL_DIM> coord,
                                                 int positionInContainer);

    /**
            Returns a boost::shared_ptr to the node object with the prescribed id.

            @param node_id integer.
     */
    boost::shared_ptr<CaVascularNetworkNode<SPATIAL_DIM> > GetNode(int node_id);

    /**
            Returns a boost::shared_ptr to the node object with the prescribed location in the spatial
            mesh. Currently, there should only ever be one node in any one location in the spatial mesh.

            To check whether there is a node at a location at all the method
            NumberOfNodesPresentAtLocation(location) may be used.
     */
    boost::shared_ptr<CaVascularNetworkNode<SPATIAL_DIM> > GetNode(ChastePoint<SPATIAL_DIM> location);

    /**
            Return mean length of vessels in network.
     */
    double GetMeanVesselLengthOfNeovasculature();

    /**
            Return mean radius of vessels in network.
     */
    double GetMeanVesselRadiusOfNeovasculature();

    /**
            Return mean tortuosity of vessels in network. Vessels whose tortuosity are infinite (i.e. self-loops)
            are ignored in this calculation.
     */
    double GetMeanVesselTortuosityOfNeovasculature();

    /**
            Return number of vessels which have length within a specified range (lowerBoundLength < length <= upperBoundLength).
     */
    int GetNumberOfVesselsByLength(double lowerBoundLength, double upperBoundLength);

    /**
            Return number of vessels which have radius within a specified range (lowerBoundRadius < length <= upperBoundRadius).
     */
    int GetNumberOfVesselsByRadius(double lowerBoundRadius, double upperBoundRadius);

    /**
            Return number of vessels which have tortuosity within a specified range (lowerBoundTortuosity  < length <= upperBoundTortuosity).
     */
    int GetNumberOfVesselsByTortuosity(double lowerBoundTortuosity, double upperBoundTortuosity);

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
    std::vector<boost::shared_ptr<CaVessel<SPATIAL_DIM> > > GetVessels();

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
    bool NodePresentAtLocation(ChastePoint<SPATIAL_DIM> location);

    /**
            Returns how many nodes are present at a location. There should only ever really be at most
            one node present at any one location.  However, the VesselNetwork may enter a
            pseudo-acceptable state within some operation implementations whereby two nodes temporarily
            occupy the same location before being merged.
     */
    int NumberOfNodesPresentAtLocation(ChastePoint<SPATIAL_DIM> location);

    /**
            Adds a vessel to the VesselNetwork. A Vessel may only be added to the network if there is
            enough free space in the spatial mesh for the vessel to be added.

            The Vessel being added should not be adjoined to any other Vessel - joining to other vessels
            in the network occurs automatically within this method. Additionally, nodes defined within
            the Vessel being added are merged with existing nodes in the vessel network. Before adding
            the Vessel, therefore, the locations of the two nodes should be defined within that Vessel.

            This method also automatically adds the Vessel to the spatial mesh.
     */
    void AddVessel(boost::shared_ptr<CaVessel<SPATIAL_DIM> > vessel);

    /**
            Checks whether the prescribed vessel is contained within the vessel network. This method is
            intended to be used mainly for debugging purposes.
     */
    bool VesselIsInNetwork(boost::shared_ptr<CaVessel<SPATIAL_DIM> > vessel);

    /**
            Checks whether the prescribed node is contained within the vessel network. This method is
            intended to be used mainly for debugging purposes.
     */
    bool NodeIsInNetwork(boost::shared_ptr<CaVascularNetworkNode<SPATIAL_DIM> > node);

    /**
            Return whether two nodes are connected.

            If the two nodes are the same this function returns true.
     */
    bool Connected(boost::shared_ptr<CaVascularNetworkNode<SPATIAL_DIM> > node1, boost::shared_ptr<CaVascularNetworkNode<SPATIAL_DIM> > node2);

    /**
            Return whether a node is connected to an input node.

            If the node is an input node this function returns true.
     */
    bool ConnectedToInputNode(boost::shared_ptr<CaVascularNetworkNode<SPATIAL_DIM> > node);

    /**
            Return whether a nodes is connected to an output node.

            If the node is an output node this function returns true.
     */
    bool ConnectedToOutputNode(boost::shared_ptr<CaVascularNetworkNode<SPATIAL_DIM> > node);

    /**
            Labels the node located at the prescribed location as an input node of the network.
            There must be a node located at the location provided.
            At present an input node must only be attached to one vessel.
     */
    void SetInputNode(ChastePoint<SPATIAL_DIM> location);

    /**
            Labels the node located at the prescribed location as an output node of the network.
            There must be a node located at the location provided.
            At present an output node must only be attached to one vessel.
     */
    void SetOutputNode(ChastePoint<SPATIAL_DIM> location);

    /**
            Save the VesselNetwork data to a file, identifiable by the prescribed string, in a format
            which may be imported to ParaView (file should be a .vtk file).
     */
    void SaveVasculatureDataToFile(string filename);

    void UpdateVascularNetwork(CaBasedCellPopulation<SPATIAL_DIM>& cell_population);

};

#endif /* CAVASCULARNETWORK_HPP_ */
