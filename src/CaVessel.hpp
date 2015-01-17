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

#ifndef CAVESSEL_HPP_
#define CAVESSEL_HPP_

#include <boost/enable_shared_from_this.hpp>
#include "SmartPointers.hpp"
#include "Exception.hpp"
#include "ChastePoint.hpp"
#include "IntraVascularChemicalCollection.hpp"
#include "Concentration.hpp"
#include "CaVascularNetworkNode.hpp"

template<unsigned DIM>
class CaVessel : public boost::enable_shared_from_this<CaVessel<DIM> >
{

private:

    /**
     *  Node at one end of the Vessel.
     */
    boost::shared_ptr<CaVascularNetworkNode<DIM> > pNode1;

    /**
     *  Node at other end of the Vessel.
     */
    boost::shared_ptr<CaVascularNetworkNode<DIM> > pNode2;

    /**
     *  Whether the Vessel has an actively migrating tip cell located at the end where node1 is
     *  located.
     */
    bool mActiveTipCellAtNode1;

    /**
     * Whether the Vessel has an actively migrating tip cell located at the end where node2 is
     * located.
     */
    bool mActiveTipCellAtNode2;

    /**
     *  Collection of VesselSegments. A vessel segment is a spatially discrete location
     *  occupied by the Vessel.
     *
     *  todo Need to think about how we introduce the following with cells instead of segments. For now use a vector of locations.
     */
    std::vector<ChastePoint<DIM> > mVesselSegmentLocations;

    /**
     * Collection of radii at vessel segment locations
     */
    std::vector<double> mVesselSegmentRadii;

    /**
     *  Time that the Vessel has experienced low wall shear stress.
     */
    double mTimeWithLowWallShearStress;

    /**
     * Radius of the vessel. Should be in units of meters.
     */
    double mRadius;

    /**
     * Previous radius of the vessel. Should be in units of meters.
     * This allows us to look for convergence of structural adaptation algorithms.
     */
    double mPreviousRadius;

    /**
     * The fractional volume of red blood cells inside a vessel.
     */
    double mHaematocritLevel;

    /**
     * The velocity of blood plasma flow through the vessel.
     */
    double mFlowVelocity;

    /**
     * The rate of blood plasma flow through the vessel.
     */
    double mFlowRate;

    /**
     * The impedance of the vessel.
     */
    double mImpedance;

    /**
     * The length of the vessel.
     */
    double mLength;

    /**
     * The shear stress experienced by the walls of the vessels.
     */
    double mWallShearStress;

    /**
     * The viscosity of the blood plasma flowing through a vessel.
     */
    double mViscosity;

    /**
     * Magnitude of the mechanical stimulus.  Many studies have shown that vessels respond
     * structurally to the mechanical forces exerted by flowing blood, i.e., transmural pressure
     * and shear stress at the endothelial surface.  Including this property allows us to model this
     * phenomenon.
     */
    double mMechanicalStimulus;

    /**
        Magnitude of the metabolic stimulus.  Inclusion of this property within a model allows us to model the fact
        that vessel segments with insufficient supply of blood may generate a signal which stimulates a
        diameter increase. If sufficiently strong, such a response has the important property of
        stabilizing  vessel network structure by preventing excessive shrinking of vessel diameters.
     */
    double mMetabolicStimulus;

    /**
        Magnitude of the shrinking stimulus.  The shrinking tendency of vessels can be interpreted
        as reflecting the basal need for factors stimulating growth of vessels to maintain or
        increase their cell mass and diameter. I.e., vessels are assumed to decrease in diameter in
        the absence of all other stimuli.
     */
    double mShrinkingStimulus;

    /**
        Magnitude of the downstream conducted stimulus. Conducted signals reaching a given vessel
        segment provide an additional stimulus to which the vessel segment may respond. This stimulus acts
        so as to increases diameter of vessel segments draining extensive networks, preventing
        low-flow and low-generation shunts.
     */
    double mDownstreamConductedStimulus;

    /**
        Magnitude of the upstream conducted stimulus. Conducted signals reaching a given vessel
        segment provide an additional stimulus to which the vessel segment may respond. This stimulus acts
        so as to increases diameter of vessel segments feeding extensive networks, preventing
        low-flow and low-generation shunts.
     */
    double mUpstreamConductedStimulus;

    /**
     * Collection of chemicals which the Vessel contains. The IntraVascularChemicalCollection
     * class allows IntraVascularChemicals to be added at will by a modeller.
     */
    IntraVascularChemicalCollection mChemicalCollection;

    /**
     *  Whether a vessel is part of the neovasculature developed since the start of a simulation.
     */
    bool mIsPartOfNeovasculature;

    /**
     * Whether the vessel is able to extend or not. This attirbute is
     * useful for implementing angiogenesis models correctly.
     */
    bool mCanExtend;

public:

    /**
       Constructor.

       When a CaVessel object is instantiated the associated CaVascularNetworkNodes have no location.
       These must be prescribed before a CaVessel is added to a CaVesselNetwork. The collection of
       segment coordinates is also empty. When segement coordinates are added to the collection of
       vessel segement coordinates collection, they must be added in order. Each coordinate must be
       located at a lattice site which neighbours the location of the previous segment coordinate
       added to the collection. A CaVessel initially contains no IntraVascularChemicals.
     */
    CaVessel();

    /**
       Destructor.
     */
    virtual ~CaVessel();

    /**
       Copies the mechanical property values and chemical concentrations from the prescribed vessel
       to this vessel. This method is useful for manipulating vessels within a CaVesselNetwork object.
       For example, this method is useful for when an existing vessel divides to form two new vessels.
       None of the spatial information about the prescribed vessel are copied (i.e. segment coordinates).
     */
    void CopyMechanicalPropertyValuesAndChemicalConcentrations(boost::shared_ptr<CaVessel<DIM> > anotherVessel);

    /**
       @return boost::shared_ptr<CaVessel<DIM> >
     */
    boost::shared_ptr<CaVessel<DIM> > shared();

    /**
       @return shared pointer to node1.
     */
    boost::shared_ptr<CaVascularNetworkNode<DIM> > GetNode1();

    /**
       @return shared pointer to node2.
     */
    boost::shared_ptr<CaVascularNetworkNode<DIM> > GetNode2();

    /**
       @return mTimeWithLowWallShearStress
     */
    double GetTimeWithLowWallShearStress();

    /**
       @return mVesselSegmentLocations
     */
    std::vector<ChastePoint<DIM> > GetSegmentCoordinates();

    /**
       @return mVesselSegmentLocations[i]
     */
    ChastePoint<DIM> GetSegmentCoordinate(unsigned i);

    /**
       @return mVesselSegmentLocations.size()
     */
    unsigned GetNumberOfSegments();

    /**
       Returns whether the Vessel has an actively migrating tip cell.

       @return true if ActiveTipCellLocatedAtNode1() or ActiveTipCellLocatedAtNode2() returns true.
     */
    bool HasActiveTipCell();

    /**
       @return mActiveTipCellAtNode1
     */
    bool ActiveTipCellLocatedAtNode1();

    /**
       @return mActiveTipCellAtNode2
     */
    bool ActiveTipCellLocatedAtNode2();

    /**
     *  @return whether vessel is attached to input node.
     */
    bool IsInputVessel();

    /**
        @return mRadius
     */
    double GetRadius();

    /**
        @return mPreviousRadius
     */
    double GetPreviousRadius();

    /**
        @return mHaematocrit
     */
    double GetHaematocritLevel();

    /**
        @return mFlowVelocity
     */
    double GetFlowVelocity();

    /**
        @return mFlowRate
     */
    double GetFlowRate();

    /**
        @return mImpedance
     */
    double GetImpedance();

    /**
        @return mLength
     */
    double GetLength();

    /**
        @return mWallShearStress
     */
    double GetWallShearStress();

    /**
        @return mViscosity
     */
    double GetViscosity();

    /**
        @return mMechanicalStimulus
     */
    double GetMechanicalStimulus();

    /**
        @return mMetabolicStimulus
     */
    double GetMetabolicStimulus();

    /**
        @return mShrinkingStimulus
     */
    double GetShrinkingStimulus();

    /**
        @return mDownstreamConductedStimulus
     */
    double GetDownstreamConductedStimulus();

    /**
        @return mUpstreamConductedStimulus
     */
    double GetUpstreamConductedStimulus();

    /*
     * todo A variant of the following should be added when cells (rather than vessel segments) are associated with vessels.
     */
//    /**
//        @return mHaematocrit
//     */
//    boost::shared_ptr<VesselSegment<DIM> > GetVesselSegment(int i);
//
//    /**
//         Returns the vessel segment at the specified index
//     */
//    boost::shared_ptr<VesselSegment<DIM> > GetVesselSegment(ChastePoint<DIM> loc);

    /**
       @return mChemicalCollection
     */
    IntraVascularChemicalCollection& rGetCollectionOfIntraVascularChemicals();

    /**
       @return the number of IntraVascularChemicals contained within the vessel.
     */
    unsigned GetNumberOfIntraVascularChemicals();

    /**
       @return the concentration of the prescribed IntraVascularChemical.
     */
    double GetIntraVascularChemicalConcentration(string chemicalname);

    /**
       @return the tortuosity of a vessel. The tortuosity is calculated as the arc-chord ratio. This is the ratio of the vessel length to the distance between the two ends of the vessel.
     */
    double GetTortuosity();

    /**
       @return mIsPartOfNeo
     */
    bool IsPartOfNeovasculature();

    /**
       @return mCanExtend
     */
    bool CanExtend();

    /**
       Assigns the concentration of the prescribed IntraVascularChemical.

       @param chemicalname
       @param concentration
     */
    void SetIntraVascularChemicalConcentration(string chemicalname, Concentration concentration);

    /**
       Assigns the prescribed node to pNode1;

       @param node node to be assigned to pNode1
     */
    void SetNode1(boost::shared_ptr<CaVascularNetworkNode<DIM> > node);

    /**
       Assigns the prescribed location to pNode1 of the vessel.

       todo Within this method the vessel is also added to the node as an adjoining vessel. Should
       look in to a way of improving how adjoining vessels are added to nodes and nodes to vessels,
       whilst maintaining a consistent state in the system.

       @param location new location of pNode1
     */
    void SetNode1Location(ChastePoint<DIM> location);

    /**
       Assigns the prescribed node to pNode2

       @param node node to be assigned to pNode2
     */
    void SetNode2(boost::shared_ptr<CaVascularNetworkNode<DIM> > node);

    /**
       Assigns the prescribed location to pNode2 of the vessel.

       todo Within this method the vessel is also added to the node as an adjoining vessel. Should
       look in to a way of improving how adjoining vessels are added to nodes and nodes to vessels,
       whilst maintaining a consistent state in the system.

       @param location new location of pNode2
     */
    void SetNode2Location(ChastePoint<DIM> location);

    /**
       Assigns a new vessel segment location to the collection of vessel segments. Vessel segments
       must be added in order and each new vessel segment must be located at a lattice site which
       neighbours the last vessel segment added to the selection.

       This method also re-calculates the length of the vessel, given that it has had
       another segment added to it.

       todo We could also potentially also automatically update the location of CaVascularNetworkNodes
       in this method. Would avoid manually updating node locations after all vessel segments have
       been added.

       @param location location of next vessel segment
     */
    void SetNextVesselSegmentCoordinate(ChastePoint<DIM> location);

    /*
     * todo Potentially add the following when cells are explicitly associated with vessels.
     */
//    /**
//            Add new VesselSegment to vessel
//     */
//    void AddVesselSegment(ChastePoint<DIM> Coords);

    /**
       Assigns the time that the Vessel has experienced low wall shear stress.

       @param value new value for mTimeWithLowWallShearStress
     */
    void SetTimeWithLowWallShearStress(double time);

    /**
       Increments the time that the Vessel has experienced low wall shear stress by the prescribed amount.

       @param value amount by which mTimeWithLowWallShearStress should be incremented
     */
    void IncrementTimeWithLowWallShearStress(double t);

    /**
       Assigns whether the Vessel has a tip cell located where node1 is located.

       @param value new value for mActiveTipCellAtNode1
     */
    void SetActiveTipCellLocatedAtNode1(bool value);

    /**
       Assigns whether the Vessel has a tip cell located where node2 is located.

       @param value new value for mActiveTipCellAtNode2
     */
    void SetActiveTipCellLocatedAtNode2(bool value);

    /**
       Assigns the prescribed value to the radius of the vessel.

       @param value new value for mRadius
     */
    void SetRadius(double value);

    /**
       Assigns the prescribed value to the previousRadius of the vessel. Note that reassignment
       of radius to previous radius must be done manually in order to avoid confusion in the
       convergence of the structural adaptation algorithm.

       @param value new value for mPreviousRadius
     */
    void SetPreviousRadius(double value);

    /**
       Assigns the prescribed value to haematocrit level.

       @param value new value for mHaematocritLevel
     */
    void SetHaematocritLevel(double value);

    /**
       Assigns the prescribed value to the velocity of blood plasma flow through the vessel.

       @param value new value for mFlowVelocity
     */
    void SetFlowVelocity(double value);

    /**
       Assigns the prescribed value to the rate of blood plasma flow through the vessel.

       @param value new value for mFlowRate
     */
    void SetFlowRate(double value);

    /**
       Assigns the prescribed value to the impedance of the vessel.

       @param value new value for mImpedance
     */
    void SetImpedance(double value);

    /**
       Assigns the prescribed value to the shear stress experienced by the walls of the vessels.

       @param value new value for mWallShearStress
     */
    void SetWallShearStress(double value);

    /**
       Assigns the prescribed value to the viscosity of the blood plasma flowing through a vessel.

       @param value new value for mViscosity
     */
    void SetViscosity(double value);

    /**
       Assigns the prescribed value to the local mechanical stimulus experienced by the vessel.

       @param value new value for mMechanicalStimulus
     */
    void SetMechanicalStimulus(double value);

    /**
       Assigns the prescribed value to the local metabolic stimulus experienced by the vessel.

       @param value new value for mMetabolicStimulus
     */
    void SetMetabolicStimulus(double value);

    /**
       Assigns the prescribed value to the shrinking stimulus experienced by the vessel.

       @param value new value for mShrinkingStimulus
     */
    void SetShrinkingStimulus(double value);

    /**
       Assigns the prescribed value to the downstream conducted stimulus experienced by the vessel.

       @param value new value for mDownstreamConductedStimulus
     */
    void SetDownstreamConductedStimulus(double value);

    /**
       Assigns the prescribed value to the upstream conducted stimulus experienced by the vessel.

       @param value new value for mUpstreamConductedStimulus
     */
    void SetUpstreamConductedStimulus(double value);

    /**
       Checks whether the prescribed node is attached to this vessel.

       @return true if vessel is attached to prescribed node
     */
    bool IsAttachedToNode(boost::shared_ptr<CaVascularNetworkNode<DIM> > node);

    /**
       Method for prescribing whether a vessel is part of the neovasculature or not.

       @param value new value for mIsPartOfNeovasculature
     */
    void SetIsPartOfNeovasculature(bool value);

    /**
       Sets whether the vessel can extend or not.

       @param value new value for mCanExtend
     */
    void CanExtend(bool value);

private:

    /**
       Assigns the prescribed value to the length of the vessel.

       @param value new value for mLength
     */
    void SetLength(double value);
};

#endif /* CAVESSEL_HPP_ */
