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

#ifndef ABSTRACTHYBRIDSOLVER_HPP_
#define ABSTRACTHYBRIDSOLVER_HPP_

#include <vector>
#include <string>
#include "OutputFileHandler.hpp"
#include "CaVascularNetwork.hpp"
#include "AbstractCellPopulation.hpp"
#include "HybridBoundaryCondition.hpp"
#include "HybridLinearEllipticPde.hpp"
#include "HybridNonLinearEllipticPde.hpp"

/**
 * An abstract solver class for hybrid continuum-discrete problems.
 * Concrete classes can solve PDEs or perform other computations based on interpolation
 * of discrete entities (points, lines) onto structured or unstructured grids.
 */
template<unsigned DIM>
class AbstractHybridSolver
{

protected:

    /**
     * The vessel network
     */
    boost::shared_ptr<CaVascularNetwork<DIM> > mpNetwork;

    /**
     * The cell population. This memory pointed to is not managed in this class.
     */
    AbstractCellPopulation<DIM>* mpCellPopulation;

    bool mCellPopulationIsSet;

    /**
     * The PDE to be solved (used in certain child classes)
     */
    boost::shared_ptr<HybridLinearEllipticPde<DIM, DIM> > mpPde;

    /**
     * The non-linear PDE to be solved (used in certain child classes)
     */
    boost::shared_ptr<HybridNonLinearEllipticPde<DIM, DIM> > mpNonLinearPde;

    /**
     * The Hybrid boundary conditions (used with PDEs in certain child classes)
     */
    std::vector<boost::shared_ptr<HybridBoundaryCondition<DIM> > > mBoundaryConditions;

    /**
     *  File handler containing output directory information
     */
    boost::shared_ptr<OutputFileHandler> mpOutputFileHandler;

    /**
     *  The filename for output
     */
    std::string mFilename;

    /**
     *  The label for the quantity being solved for
     */
    std::string mLabel;

    std::vector<double> mPointSolution;

    bool IsSetupForSolve;

    bool mWriteSolution;

public:

    /**
     * Constructor
     */
    AbstractHybridSolver();

    /**
     * Destructor
     */
    virtual ~AbstractHybridSolver();

    /**
     * Add a hybrid boundary condition for the domain
     * @param pBoundaryCondition the boundary condition
     */
    void AddBoundaryCondition(boost::shared_ptr<HybridBoundaryCondition<DIM> > pBoundaryCondition);

    void SetWriteSolution(bool write=true);

    /**
     * Return the PDE
     * @return the hybrid linear elliptic pde
     */
    boost::shared_ptr<HybridLinearEllipticPde<DIM, DIM> > GetPde();

    boost::shared_ptr<HybridNonLinearEllipticPde<DIM, DIM> > GetNonLinearPde()
    {
        return mpNonLinearPde;
    }

    /**
     * Return the solver output sample at discrete points. Different sampling strategies are implemented in child classes.
     * @param samplePoints the points to be sampled at
     * @param samplingStrategy use the default sampling strategy
     * @return a vector of the point values
     */
    virtual std::vector<double> GetSolutionAtPoints(std::vector<c_vector<double, DIM> > samplePoints, bool samplingStrategy = true) = 0;

    /**
     * Set the cell population
     * @param rCellPopulation a reference to the cell population
     */
    void SetCellPopulation(AbstractCellPopulation<DIM>& rCellPopulation);

    /**
     *  Set the PDE to be solved in certain child classes
     * @param pPde the pde to be solved
     */
    void SetPde(boost::shared_ptr<HybridLinearEllipticPde<DIM, DIM> > pPde);

    void SetNonLinearPde(boost::shared_ptr<HybridNonLinearEllipticPde<DIM, DIM> > pPde)
    {
        mpNonLinearPde = pPde;
    }

    void SetLabel(const std::string& label);

    const std::string& GetLabel();

    /**
     * Operations to be performed prior to the first solve
     */
    virtual void Setup() = 0;

    /**
     * Set the vessel network
     * @param pNetwork the vessel network
     */
    void SetVesselNetwork(boost::shared_ptr<CaVascularNetwork<DIM> > pNetwork);

    /**
     * Do the solve
     * @param writeSolution whether to write the solution to file
     */
    virtual void Solve() = 0;

    /**
     * Set the file handler containing the working directory
     * @param pOutputFileHandler the file handler containing the working directory
     */
    void SetFileHandler(boost::shared_ptr<OutputFileHandler> pOutputFileHandler);

    /**
     * Set the file name for output
     * @param filename the file name
     */
    void SetFileName(const std::string& rFilename);

    virtual bool HasRegularGrid();

    /**
     * Operations to be performed prior to every solve
     */
    virtual void Update() = 0;

    virtual void UpdateCellData() = 0;

    std::vector<double> GetPointSolution();
    /**
     * Write the solution to file
     */
    virtual void Write() = 0;
};

#endif /* ABSTRACTHYBRIDSOLVER_HPP_ */
