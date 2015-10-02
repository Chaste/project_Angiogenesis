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

#ifndef SIMPLEFLOWSOLVER_HPP_
#define SIMPLEFLOWSOLVER_HPP_

#include <boost/shared_ptr.hpp>

#include "CaVascularNetwork.hpp"
#include "CaVessel.hpp"
#include "VascularNode.hpp"
#include "LinearSystem.hpp"

template<unsigned DIM>
class SimpleFlowSolver
{
private:

    std::vector<boost::shared_ptr<VascularNode<DIM> > > mNodes;
    std::vector<boost::shared_ptr<CaVessel<DIM> > > mVessels;
    boost::shared_ptr<CaVascularNetwork<DIM> > mpVesselNetwork;
    std::vector<std::vector<unsigned> > mNodeVesselConnectivity;
    std::vector<std::vector<unsigned> > mNodeNodeConnectivity;
    std::vector<unsigned> mBoundaryConditionNodeIndices;
    std::vector<unsigned> mUnconnectedNodeIndices;
    boost::shared_ptr<LinearSystem> mpLinearSystem;
    bool mUseDirectSolver;
    double mMultiplier;
    bool mIsSetUp;

public:

    /**
     * Constructor.
     */
    SimpleFlowSolver();

    /*
     * Construct a new instance of the class and return a shared pointer to it.
     */
    static boost::shared_ptr<SimpleFlowSolver<DIM> > Create();

    void SetUseDirectSolver(bool useDirectSolver);

    /**
     * Destructor.
     */
    ~SimpleFlowSolver();

    bool IsSetUp();

    /**
     * Set up the flow solver;
     */
    void SetUp(boost::shared_ptr<CaVascularNetwork<DIM> > pVascularNetwork =
            boost::shared_ptr<CaVascularNetwork<DIM> >());

    /**
     * Update the impedances in the system matrix
     */
    void UpdateImpedances();

    void SetImpedanceScaleFactor(double multiplier);

    /**
     * Implement flow solver;
     */
    void Implement(
            boost::shared_ptr<CaVascularNetwork<DIM> > pVascularNetwork = boost::shared_ptr<CaVascularNetwork<DIM> >());

};

#endif /* SIMPLEFLOWSOLVER_HPP_ */
