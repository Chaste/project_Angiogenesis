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

#ifndef VASCULARTUMOURSOLVER_HPP_
#define VASCULARTUMOURSOLVER_HPP_

#include <vector>
#include <string>

#include "StructuralAdaptationSolver.hpp"
#include "SmartPointers.hpp"
#include "CaVascularNetwork.hpp"
#include "AbstractHybridSolver.hpp"
#include "FlowSolver.hpp"
#include "AbstractCellPopulation.hpp"
#include "AngiogenesisSolver.hpp"

/**
 * This class manages the solution of vascular tumour growth problems. It steps through time,
 * solves a collection of hybrid discrete-continuum PDEs and, if required, updates the vessel network. For linking
 * with discrete cell models it can be added to a VascularTumourGrowth simulation modifier.
 */
template<unsigned DIM>
class VascularTumourSolver
{
    /**
     * The vessel network
     */
    boost::shared_ptr<CaVascularNetwork<DIM> > mpNetwork;

    /**
     * The end time for solves if used in standalone
     */
    double mEndTime;

    /**
     * The frequency of file based output
     */
    unsigned mOutputFrequency;

    /**
     * Filehandler containing output directory information
     */
    boost::shared_ptr<OutputFileHandler> mpOutputFileHandler;

    /**
     * The collection of hybrid solvers
     */
    std::vector<boost::shared_ptr<AbstractHybridSolver<DIM> > > mHybridSolvers;

    /**
     * The flow solver in the vessel network
     */
    boost::shared_ptr<FlowSolver<DIM> > mpFlowSolver;

    /**
     * The structural adaptation solver for the vessel network
     */
    boost::shared_ptr<StructuralAdaptationSolver<DIM> > mpStructuralAdaptationSolver;

    /**
     * The angiogenesis solver for the vessel network
     */
    boost::shared_ptr<AngiogenesisSolver<DIM> > mpAngiogenesisSolver;

public:

    /**
     * Constructor.
     */
    VascularTumourSolver();

    /**
     * Destructor.
     */
    virtual ~VascularTumourSolver();

    /**
     * Factory constructor method
     * @return a shared pointer to a new solver
     */
    static boost::shared_ptr<VascularTumourSolver> Create();

    /**
     * Add a hybrid solver to the collection
     * @param pHybridSolver a hybrid discrete-continuum solver
     */
    void AddHybridSolver(boost::shared_ptr<AbstractHybridSolver<DIM> > pHybridSolver);

    /**
     * Return the current hybrid solvers
     * @return the hybrid solvers
     */
    std::vector<boost::shared_ptr<AbstractHybridSolver<DIM> > > GetHybridSolvers();

    /**
     * Set the angiogenesis solver for the network
     * @param pAngiogenesisSolver the solver for structural adaptation
     */
    void SetAngiogenesisSolver(boost::shared_ptr<AngiogenesisSolver<DIM> > pAngiogenesisSolver);

    /**
     * Set the simulation end time
     * @param time the simulation end time
     */
    void SetEndTime(double time);

    /**
     * Set the flow solver for the network
     * @pFlowSolver the flow solver for the network
     */
    void SetFlowSolver(boost::shared_ptr<FlowSolver<DIM> > pFlowSolver);

    /**
     * Set the output directory for results
     * @param pFileHandler output file handler containing output directory information
     */
    void SetOutputFileHandler(boost::shared_ptr<OutputFileHandler> pFileHandler);

    /**
     * Set the results output frequency
     * @param frequency the frequency of simulaiton output
     */
    void SetOutputFrequency(unsigned frequency);

    /**
     * Set the structural adaptation solver for the network
     * @param pStructuralAdaptationSolver the solver for structural adaptation
     */
    void SetStructuralAdaptationSolver(boost::shared_ptr<StructuralAdaptationSolver<DIM> > pStructuralAdaptationSolver);

    void SetupFromModifier(AbstractCellPopulation<DIM,DIM>& rCellPopulation, const std::string& rDirectory);

    void Setup();

    void UpdateCellData(std::vector<std::string> labels);

    /**
     * Set the vessel network
     * @param pNetwork the vessel network
     */
    void SetVesselNetwork(boost::shared_ptr<CaVascularNetwork<DIM> > pNetwork);

    /**
     * Increment one step in time
     */
    void Increment();

    /**
     * Run until the specified end time
     */
    void Run();
};

#endif /* VASCULARTUMOURSOLVER_HPP_ */