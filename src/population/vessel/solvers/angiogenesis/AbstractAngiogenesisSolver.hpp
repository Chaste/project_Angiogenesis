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

#ifndef ABSTRACTANGIOGENESISSOLVER_HPP_
#define ABSTRACTANGIOGENESISSOLVER_HPP_

#include <vector>
#include <string>

#include "CaVascularNetwork.hpp"
#include "AbstractHybridSolver.hpp"
#include "SmartPointers.hpp"

template<unsigned DIM>
class AbstractAngiogenesisSolver
{

    boost::shared_ptr<CaVascularNetwork<DIM> > mpNetwork;
    double mGrowthVelocity;
    double mTimeIncrement;
    double mEndTime;
    unsigned mOutputFrequency;
    std::string mOutputDirectory;
    double mNodeAnastamosisRadius;
    std::vector<boost::shared_ptr<AbstractHybridSolver<DIM> > > mPdeSolvers;
    bool mSolveFlow;
    double mSproutingProbability;

public:

    /**
     * Constructor.
     * @param pNetwork the network to perform angiogenesis on
     */
    AbstractAngiogenesisSolver(boost::shared_ptr<CaVascularNetwork<DIM> > pNetwork);

    /**
     * Destructor.
     */
    virtual ~AbstractAngiogenesisSolver();

    virtual c_vector<double, DIM> GetGrowthDirection(c_vector<double, DIM> currentDirection);

    void AddPdeSolver(boost::shared_ptr<AbstractHybridSolver<DIM> > pPdeSolver);

    void SetSolveFlow(bool solveFlow=true);

    void SetSproutingProbability(double sproutingProbability);

    void SetOutputDirectory(const std::string& rDirectory);

    std::vector<boost::shared_ptr<AbstractHybridSolver<DIM> > > GetPdeSolvers();

    void UpdateNodalPositions(const std::string& speciesLabel = "Default");

    void DoSprouting();

    void DoAnastamosis();

    void Increment();

    /**
     * Run the solver.
     */
    void Run();
};

#endif /* ABSTRACTANGIOGENESISSOLVER_HPP_ */
