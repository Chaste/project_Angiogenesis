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

#ifndef SIMPLECELLPOPULATION_HPP_
#define SIMPLECELLPOPULATION_HPP_

#include <vector>
#define _BACKWARD_BACKWARD_WARNING_H 1 //Cut out the vtk deprecated warning
#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include "Part.hpp"
#include "SmartPointers.hpp"
#include "CaVascularNetwork.hpp"
#include "SimpleCell.hpp"

/* A minimal cell population class.
 */

template<unsigned DIM>
class SimpleCellPopulation
{
    std::vector<boost::shared_ptr<SimpleCell<DIM> > > mCells;

public:

    /* Constructor
     */

    SimpleCellPopulation();

    /* Factory constructor method
     * @return a shared pointer to a new population
     */
    static boost::shared_ptr<SimpleCellPopulation<DIM> > Create();

    /* Desctructor
     */
    ~SimpleCellPopulation();

    std::vector<boost::shared_ptr<SimpleCell<DIM> > > GetCells();

    void AddCell(boost::shared_ptr<SimpleCell<DIM> > pCell);

    void AddCells(std::vector<boost::shared_ptr<SimpleCell<DIM> > > cells);

    void BooleanWithVesselNetwork(boost::shared_ptr<CaVascularNetwork<DIM> > pNetwork);

    void GenerateCellsOnGrid(unsigned xDim = 10, unsigned yDim = 10, unsigned zDim = 10,
                             double spacing =1.0, c_vector<double, DIM> origin = zero_vector<double>(DIM));

    void GenerateCellsOnGrid(boost::shared_ptr<Part<DIM> > pPart, double spacing =1.0);

    vtkSmartPointer<vtkPoints> GetVtk();

    void Write(const std::string& rFileName);
};

#endif /* SIMPLECELLPOPULATION_HPP_*/