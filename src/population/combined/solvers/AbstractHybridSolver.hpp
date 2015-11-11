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
#define _BACKWARD_BACKWARD_WARNING_H 1 //Cut out the vtk deprecated warning for now (gcc4.3)
#include <vtkImageData.h>
#include <vtkSmartPointer.h>
#include "HybridLinearEllipticPde.hpp"
#include "DirichletBoundaryCondition.hpp"

/*
 * An abstract solver class for linear elliptic PDEs which can include
 * discrete representations of cells and vessels.
 *
 */
template<unsigned DIM>
class AbstractHybridSolver
{

protected:

    /* The vessel network
    */
    boost::shared_ptr<CaVascularNetwork<DIM> > mpNetwork;

    /* The pde to be solved
    */
    boost::shared_ptr<HybridLinearEllipticPde<DIM, DIM> > mpPde;

    /* The dirichlet boundary conditions
    */
    std::vector<boost::shared_ptr<DirichletBoundaryCondition<DIM> > > mDirichletBoundaryConditions;

    /* Filehandler containing output directory info
    */
    boost::shared_ptr<OutputFileHandler> mpOutputFileHandler;

    /* The filename for output
    */
    std::string mFilename;

public:

    /* Constructor
     */
    AbstractHybridSolver();

    /* Destructor
     */
    virtual ~AbstractHybridSolver();

    /* Get the solution as vtk image data. If the solution is on an unstructured grid it must be first resampled.
     * @return the solution as vtk image data
     */
    virtual vtkSmartPointer<vtkImageData> GetVtkSolution();

    boost::shared_ptr<HybridLinearEllipticPde<DIM, DIM> > GetPde();

    virtual std::vector<double> GetSolutionAtPoints(std::vector<c_vector<double, DIM> > samplePoints, bool useVtkSampling = true);

    /* Set a boundary condition for the domain
     * @param pBoundaryCondition the boundary condition
     */
    void AddDirichletBoundaryCondition(boost::shared_ptr<DirichletBoundaryCondition<DIM> > pBoundaryCondition);

    /* Set the pde to be solved
     * @param pPde the pde to be solved
     */
    void SetPde(boost::shared_ptr<HybridLinearEllipticPde<DIM, DIM> > pPde);

    /* Solve the pde
     * @param writeSolution whether to write the solution to file
     */
    virtual void Solve(bool writeSolution = false) = 0;

    virtual void Setup() = 0;

    /* Set the working directory
     * @param directory the working directory
     */
    void SetFileHandler(boost::shared_ptr<OutputFileHandler> pOutputFileHandler);

    /* Set the file name for output
     * @param filename the file name
     */
    void SetFileName(const std::string& filename);

protected:

    /* Write the solution to file
     */
    void Write();

    /* Update the solution
     * @param data solution data map
     */
    void UpdateSolution(std::map<std::string, std::vector<double> >& data);
};

#endif /* ABSTRACTHYBRIDSOLVER_HPP_ */
