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

#ifdef CHASTE_ANGIOGENESIS_PYTHON
#ifndef PDEWRAPPERS_HPP_
#define PDEWRAPPERS_HPP_

#include <boost/python.hpp>
#include <boost/python/module.hpp>
#include <vtkSmartPointer.h>
#include <vtkImageData.h>
#include "SmartPointers.hpp"
#include "FunctionMap.hpp"
#include "AbstractRegularGridHybridSolver.hpp"
#include "AbstractHybridSolver.hpp"
#include "RegularGrid.hpp"

using namespace boost::python;

typedef AbstractHybridSolver<3> AbstractHybridSolver3;
struct AbstractHybridSolverWrap : AbstractHybridSolver3, wrapper<AbstractHybridSolver3>
{
    void Solve()
    {
        this->get_override("Solve")();
    }

    void Update()
    {
        this->get_override("Update")();
    }

    void Setup()
    {
        this->get_override("Setup")();
    }

    void Write()
    {
        this->get_override("Write")();
    }

    std::vector<double> GetSolutionAtPoints(std::vector<c_vector<double, 3> > samplePoints)
    {
        this->get_override("GetSolutionAtPoints")();
    }

    std::vector<double> GetSolutionAtGridPoints(boost::shared_ptr<RegularGrid<3> > pGrid)
    {
        this->get_override("GetSolutionAtGridPoints")();
    }

    void UpdateCellData()
    {
        this->get_override("UpdateCellData")();
    }
};

typedef AbstractRegularGridHybridSolver<3> AbstractRegularGridHybridSolver3;
struct AbstractRegularGridHybridSolverWrap : AbstractRegularGridHybridSolver3, wrapper<AbstractRegularGridHybridSolver3>
{
    void Solve()
    {
        this->get_override("Solve")();
    }

    void Update()
    {
        this->get_override("Update")();
    }

    vtkSmartPointer<vtkImageData> GetVtkSolution()
    {
        if (override GetVtkSolution = this->get_override("GetVtkSolution"))
        {
            return GetVtkSolution();
        }
        return AbstractRegularGridHybridSolver<3>::GetVtkSolution();
    }

    vtkSmartPointer<vtkImageData> default_GetVtkSolution()
    {
        return this->AbstractRegularGridHybridSolver<3>::GetVtkSolution();
    }
};

#endif /* PDEWRAPPERS_HPP_ */
#endif // CHASTE_ANGIOGENESIS_PYTHON
