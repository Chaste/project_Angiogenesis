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

#ifndef ABSTRACTREGULARGRIDHYBRIDSOLVER_HPP_
#define ABSTRACTREGULARGRIDHYBRIDSOLVER_HPP_

#include <vector>
#include <string>
#define _BACKWARD_BACKWARD_WARNING_H 1 //Cut out the vtk deprecated warning for now (gcc4.3)
#include <vtkImageData.h>
#include <vtkSmartPointer.h>
#include "Part.hpp"
#include "AbstractHybridSolver.hpp"
#include "RegularGrid.hpp"

/*
 * An abstract solver class for linear elliptic PDEs on regular grids which can include
 * discrete representations of cells and vessels.
 */
template<unsigned DIM>
class AbstractRegularGridHybridSolver : public AbstractHybridSolver<DIM>
{

protected:

    /* The solution in the form of vtk image data
    */
    vtkSmartPointer<vtkImageData> mpRegularGridVtkSolution;

    /**
     * The grid for solvers using regular grids
     */
    boost::shared_ptr<RegularGrid<DIM> > mpRegularGrid;

public:

    /* Constructor
     */
    AbstractRegularGridHybridSolver();

    /* Destructor
     */
    virtual ~AbstractRegularGridHybridSolver();

    boost::shared_ptr<RegularGrid<DIM> > GetGrid();

    std::vector<double> GetSolutionAtPoints(std::vector<c_vector<double, DIM> > samplePoints, bool useVtkSampling = true);

    std::vector<double> GetSolutionOnLine(double sampleSpacing,
                                          c_vector<double, DIM> start_point,
                                          c_vector<double, DIM> end_point,
                                          const std::string& rSpeciesLabel = "Default");

    vtkSmartPointer<vtkImageData> GetSolutionOnVolume(double sampleSpacing = 1.0,
                                            std::vector<unsigned> dimensions = std::vector<unsigned>(DIM,10),
                                            const std::string& rSpeciesLabel = "Default",
                                            c_vector<double, DIM> origin = zero_vector<double>(DIM));

    /* Get the solution as vtk image data
     * @return the solution as vtk image data
     */
    vtkSmartPointer<vtkImageData> GetVtkSolution();

    double GetVolumeAverageSolution(const std::string& arrayName);

    /* Set the number of grid points in each dimension based on the bounding box of
     * a part.
     * @param pPart the part from which to get the dimension
     */
    void SetGridFromPart(boost::shared_ptr<Part<DIM> > pPart, double gridSize);

    void SetGrid(boost::shared_ptr<RegularGrid<DIM> > pRegularGrid);

    virtual void Setup();

    virtual void Solve(bool writeSolution = false) = 0;

    void WriteHistograms(const std::string& arrayName, const std::string& fileName, double binSize, unsigned numberOfBins);

protected:

    /* Write the solution to file
     */
    void Write();

    /* Update the solution
     * @param data solution data map
     */
    void UpdateSolution(std::map<std::string, std::vector<double> >& data);
};

#endif /* ABSTRACTREGULARGRIDHYBRIDSOLVER_HPP_ */
