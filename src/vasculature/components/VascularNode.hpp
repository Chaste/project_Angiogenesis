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
#include <boost/weak_ptr.hpp>
#include <boost/enable_shared_from_this.hpp>
#include "ChastePoint.hpp"
#include "Exception.hpp"
#include "Cell.hpp"
#include "AbstractCellPopulation.hpp"
#include "SmartVasculaturePointers.hpp"
#include "VasculatureData.hpp"

/**
 *  Forward declaration to allow segments to manage adding and removing themselves from nodes.
 */
template<unsigned DIM>
class CaVesselSegment;

/*
 * Nodes are point locations along a vessel. They are useful for describing the positions of
 * vessel segments.
 */

template<unsigned DIM>
class VascularNode : public boost::enable_shared_from_this<VascularNode<DIM> >
{

    /**
     *  Allow segments to manage adding and removing themselves from nodes.
     */
    friend class CaVesselSegment<DIM>;

private:

    /**
     *   Location of a node. Used if a CellPtr has not been assigned.
     */
    ChastePoint<DIM> mLocation;

    /**
     *   Pointer to an associated Cell.
     */
    CellPtr mpCell;

    /**
     *   Pointer to an associated Cell Population.
     */
    boost::shared_ptr<AbstractCellPopulation<DIM> > mpCellPopulation;

    /**
     * Container for non-spatial node data.
     */
    VasculatureData mDataContainer;

    /**
     * Id tag, can be useful for storing segment-node relationships in the VesselNetwork class.
     */
    unsigned mId;

    /**
     * Label tag, can be useful for identifying input and output nodes.
     */
    std::string mLabel;

    /**
     *   Collection of pointers to Vessel Segments connected to this node.
     */
    std::vector<boost::weak_ptr<CaVesselSegment<DIM> > > mVesselSegments;

public:

    /*
     * Constructor
     *
     * A node must be specified with a location.
     */
    VascularNode(const ChastePoint<DIM>& rLocation);

    /*
     * Constructor
     *
     * A node must be specified with a location.
     */
    VascularNode(double point1, double point2, double point3 = 0.0);
    
    /*
     * Destructor
     */
    ~VascularNode();

    /*
     * Construct a new instance of the class and return a shared pointer to it.
     */
    static boost::shared_ptr<VascularNode<DIM> > Create(const ChastePoint<DIM>& rLocation);

    /*
     * Construct a new instance of the class and return a shared pointer to it.
     */
    static boost::shared_ptr<VascularNode<DIM> > Create(double point1, double point2, double point3 = 0.0);

    /**
     *  Return a pointer to the associated Cell.
     *
     *  @return mpCell
     */
    CellPtr GetCell() const;
    
    /**
     *  Return the node data for the input key. An attempt is made
     *  to cast to type T.
     *  @return T data
     */
    template<typename T> T GetData(const std::string& rKey);

    /**
     *  Return a const reference to the non-spatial data container.
     *
     *  @return mDataContainer
     */
    const VasculatureData& rGetDataContainer();   
    
    /**
     *  Return a vector of data keys for the node. Input true if
     *  the corresponding value should be castable to double.
     *
     *  @return std::vector<std::string>
     */
    std::vector<std::string> GetDataKeys(bool castable_to_double = false) const;   
    
    /**
     *  Return the distance between the input point and the node
     *
     *  @return double
    */
    double GetDistance(const ChastePoint<DIM>& rPoint);

    /**
     *  Return the Id
     *
     *  @return mId
     */
    unsigned GetId() const;

    /**
     *  Return a const reference to the Label
     *
     *  @return mLabel
     */
    const std::string& rGetLabel() const;

    /**
     *  Return the location of the node or, if there is one, the associated Cell.
     *
     *  @return mLocation
     */
    ChastePoint<DIM> GetLocation() const;

    /**
     *  Return the number of attached segments
     */
    unsigned GetNumberOfSegments() const;

    /**
     *  Return a boost::shared_ptr to VesselSegment index.
     */
    boost::shared_ptr<CaVesselSegment<DIM> > GetVesselSegment(unsigned index) const;

    /**
     *  Return a vector of boost::shared_ptr to the VesselSegments.
     */
    std::vector<boost::shared_ptr<CaVesselSegment<DIM> > > GetVesselSegments() const;

    /**
     *  Return true if there is an associated Cell.
     *
     *  @return bool
     */
    bool HasCell() const;
    
     /**
     *  Return true if the node has data corresponding to the input key.
     *
     *  @return bool
     */
    bool HasDataKey(const std::string& rKey) const;

    /**
     *  Return true if the input segment is attached to the node
     *
     *  @return bool
     */
    bool IsAttachedTo(const boost::shared_ptr<CaVesselSegment<DIM> > pSegment) const;

    /**
     *  Return true if the node is coincident with the input location
     *
     *  @return bool
     */
    bool IsCoincident(const ChastePoint<DIM>& rPoint) const;

    /**
     *  Return true if the node is coincident with the input node
     *
     *  @return bool
     */
    bool IsCoincident(const boost::shared_ptr<VascularNode<DIM> > node) const;

    /**
     *  Assign a Cell to the node. Throw an Exception if the Cell is not a member of the
     *  Cell Population or if a Cell Population has not been set. Overwrite any existing Cell.
     */
    void SetCell(CellPtr pCell);

    /**
     *  Assign a Cell Population to the node. If an existing Cell is not a member of the population
     *  remove it.
     */
    void SetCellPopulation(boost::shared_ptr<AbstractCellPopulation<DIM> > pCellPopulation);

    /**
     *  Add data of any type to the node using the identifying key
     */
    template<typename T> void SetData(const std::string& rKey, T value);

    /**
     *  Over-write the node's non-spatial DataContainer
     *
     *  This can be useful when copying data from an existing node.
     */
    void SetDataContainer(const VasculatureData& rDataContainer);

    /**
     *  Assign the Id
     *
     */
    void SetId(unsigned id);

    /**
     *  Assign the Label
     *
     */
    void SetLabel(const std::string& rLabel);

    /**
     *  Set the location of the node. This breaks any links with an assigned Cell, so if there is an
     *  assigned Cell remove it.
     */
    void SetLocation(const ChastePoint<DIM>& rLocation);

    /**
     *  Remove the assigned Cell
     */
    void RemoveCell();

private:

    /**
       Adds an adjoining Vessel segment the node.
     */
    void AddSegment(boost::shared_ptr<CaVesselSegment<DIM> > pVesselSegment);

    /**
       Removes an adjoining vessel segment from to node.
     */
    void RemoveSegment(boost::shared_ptr<CaVesselSegment<DIM> > pVesselSegment);

};

#endif /* VASCULARNODE_HPP_ */
