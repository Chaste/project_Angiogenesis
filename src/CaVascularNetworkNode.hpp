/*
 * CaVascularNetworkNode.hpp
 *
 *  Created on: 13 Jan 2015
 *      Author: connor
 */

#ifndef CAVASCULARNETWORKNODE_HPP_
#define CAVASCULARNETWORKNODE_HPP_

#include <ChastePoint.hpp>

#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>
#include <boost/enable_shared_from_this.hpp>

#include <vector>

// forward declaration of vessel class - vascular network node referenced in vessel class
template <unsigned SPATIAL_DIM>
class CaVessel;

template<unsigned SPATIAL_DIM>
class CaVascularNetworkNode : public boost::enable_shared_from_this<CaVascularNetworkNode<SPATIAL_DIM> >
{
private:

    /**
        Location of node.
     */
    ChastePoint<SPATIAL_DIM> mLocation;


    /**
        Pressure in VesselNetwork at node location. Flow through Vessels in a VesselNetwork is dictated by
        the pressure differences between nodes in the network.

        Pressure should be in units of Pa.
     */
    double mPressure;

    /**
        Collection of CaVessel objects which are adjoint to this node.
     */
    std::vector<boost::weak_ptr<CaVessel<SPATIAL_DIM> > > pAdjoiningVessels;

    /**
        Whether the node is an input node to a CaVesselNetwork or not.
     */
    bool mIsInputNode;

    /**
        Whether the node is an output node of a CaVesselNetwork or not.
     */
    bool mIsOutputNode;

public:

    /*
     * Constructor
     *
     * Upon instantiation the node has no prescribed location, no adjoining vessels, pressure at the
     * node is zero and the node is neither an input or output node to a network.
     *
     */
    CaVascularNetworkNode();

    /*
     * Destructor
     */
    virtual ~CaVascularNetworkNode();

    /**
       Returns a boost::shared_ptr to this CaVascularNetworkNode.

       @return boost::shared_ptr<CaVascularNetworkNode<DIM,T> >
     */
    boost::shared_ptr<CaVascularNetworkNode<SPATIAL_DIM> > shared();

    /**
       Returns the location of the node.

       @return mLocation
     */
    ChastePoint<SPATIAL_DIM> GetLocation();

    /**
       Returns the numerical value of the pressure at the node.

       @return mPressure
     */
    double GetPressure();

    /**
       Returns the number of vessels which are adjoint to the node.

       @return pAdjoiningVessels.size()
     */
    unsigned GetNumberOfAdjoiningVessels();

    /**
       Returns a boost::shared_ptr to Vessel i which is adjoint to this node.

       @return pAdjoiningVessels[i]
     */
    boost::shared_ptr<CaVessel<SPATIAL_DIM> > GetAdjoiningVessel(unsigned i);

    /**
       Returns whether the node is an input node to a VesselNetwork.

       @return mIsInputNode
     */
    bool IsInputNode();

    /**
       Returns whether the node is an output node to a VesselNetwork.

       @return mIsInputNode
     */
    bool IsOutputNode();

    /**
       Assigns the prescribed location to the node.

       @param loc new location of node
     */
    void SetLocation(ChastePoint<SPATIAL_DIM> loc);


    /**
       Assigns the prescribed pressure to the node.

       @param pressure new pressure at node
     */
    void SetPressure(double pressure);

    /**
       Adds an adjoining Vessel to the node.
       The prescribed Vessel may be attached to the same node twice if it loops around on itself.
       In this case the node must be both node1 and node2 of the Vessel.

       todo: at some point should possibly check that the vessel being adjoined to node has not got
       an actively migrating tip located at this node.

       @param vessel vessel to be added to node
     */
    void AddAdjoiningVessel(boost::shared_ptr<CaVessel<SPATIAL_DIM> > vessel);

    /**
       Removes an adjoining vessel from to node.

       @param vessel vessel to be removed from node
     */
    void RemoveAdjoiningVessel(boost::shared_ptr<CaVessel<SPATIAL_DIM> > vessel);

    /**
       Checks whether the prescribed vessel is attached to this node.

       @param vessel vessel which may or may not be attached to node
     */
    bool IsAttachedToVessel(boost::shared_ptr<CaVessel<SPATIAL_DIM> > vessel);

    /**
       Assigns whether the node is an input node to a vessel network.

       @param value whether node is input node
     */
    void SetIsInputNode(bool value);

    /**
       Assigns whether the node is an output node to a vessel network.

       @param value whether node is input node
     */
    void SetIsOutputNode(bool value);

};

#endif /* CAVASCULARNETWORKNODE_HPP_ */
