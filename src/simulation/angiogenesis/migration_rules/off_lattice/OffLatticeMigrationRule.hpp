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

#ifndef OffLatticeMigrationRule_HPP_
#define OffLatticeMigrationRule_HPP_

#include <vector>
#include <string>
#include "AbstractMigrationRule.hpp"
#include "VesselNode.hpp"
#include "SmartPointers.hpp"

/**
 * An off-lattice migration rule for tip cells
 */
template<unsigned DIM>
class OffLatticeMigrationRule : public AbstractMigrationRule<DIM>
{
    /**
     * Global direction vectors, x (1,0,0)
     */
    c_vector<double, DIM> mGlobalX;

    /**
     * Global direction vectors, y (0,1,0)
     */
    c_vector<double, DIM> mGlobalY;

    /**
     * Global direction vectors, z (0,0,1)
     */
    c_vector<double, DIM> mGlobalZ;

    /**
     * Mean angle between current and new directions about global axes
     */
    std::vector<double> mMeanAngles;

    /**
     * Deviation in angle between current and new directions about global axes
     */
    std::vector<double> mSdvAngles;

    /**
     * Tip cell velocity
     */
    double mVelocity;

    double mChemotacticStrength;

    double mAttractionStrength;

    /**
     * Length of probe into solution
     */
    double mProbeLength;


    bool mIsSprouting;

public:

    /**
     * Constructor.
     */
    OffLatticeMigrationRule();

    /**
     * Destructor.
     */
    virtual ~OffLatticeMigrationRule();

    /*
     * Construct a new instance of the class and return a shared pointer to it.
     * @return pointer to a new class instance
     */
    static boost::shared_ptr<OffLatticeMigrationRule<DIM> > Create();

    void SetIsSprouting(bool isSprouting = true);

    void SetSproutingVelocity(double velocity);

    void SetChemotacticStrength(double strength);

    void SetAttractionStrength(double strength);

    /**
     * Return the movement vector (new_location - oriringal_location) for the input nodes, if they can't move set it to the zero vector
     * @param rNodes nodes to calculate indices
     * @return a vector of movement vectors
     */
    std::vector<c_vector<double, DIM> > GetDirections(const std::vector<boost::shared_ptr<VesselNode<DIM> > >& rNodes);


    std::vector<c_vector<double, DIM> > GetDirectionsForSprouts(const std::vector<boost::shared_ptr<VesselNode<DIM> > >& rNodes);
};

#endif /* OffLatticeMigrationRule_HPP_ */
