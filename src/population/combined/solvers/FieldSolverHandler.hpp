///*
//
//Copyright (c) 2005-2015, University of Oxford.
//All rights reserved.
//
//University of Oxford means the Chancellor, Masters and Scholars of the
//University of Oxford, having an administrative office at Wellington
//Square, Oxford OX1 2JD, UK.
//
//This file is part of Chaste.
//
//Redistribution and use in source and binary forms, with or without
//modification, are permitted provided that the following conditions are met:
// * Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
// * Neither the name of the University of Oxford nor the names of its
//   contributors may be used to endorse or promote products derived from this
//   software without specific prior written permission.
//
//THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
//AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
//IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
//ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
//LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
//CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
//GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
//HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
//LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
//OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
//*/
//
//#ifndef FIELDSOLVEREHANDLER_HPP_
//#define FIELDSOLVEREHANDLER_HPP_
//
//#include <map>
//#include <memory>
//
//#include "ChasteSerialization.hpp"
//#include <boost/serialization/vector.hpp>
//
//#include "AbstractCellPopulation.hpp"
//#include "PdeAndBoundaryConditions.hpp"
//#include "BoundaryConditionsContainer.hpp"
//#include "TetrahedralMesh.hpp"
//#include "ChasteCuboid.hpp"
//
///**
// * Class for handling the solution of field problems with cells and vessels.
// * Can include solving PDEs or also more generic field calculations, such as distance maps.
// */
//template<unsigned DIM>
//class FieldSolverHandler : public CellBasedPdeHandler
//{
//
//public:
//
//    /**
//     * Constructor.
//     *
//     * @param pCellPopulation pointer to a cell population
//     * @param deleteMemberPointersInDestructor whether to delete member pointers in the destructor (defaults to false)
//     */
//    FieldSolverHandler(AbstractCellPopulation<DIM>* pCellPopulation, bool deleteMemberPointersInDestructor=false);
//
//    /**
//     * Destructor.
//     */
//    virtual ~FieldSolverHandler();
//
//    /**
//     * Solve the PDE and write the solution to file.
//     *
//     * @param samplingTimestepMultiple the ratio of the number of actual timesteps to the number of timesteps
//     *     at which results are written to file.
//     */
//    virtual void SolvePdeAndWriteResultsToFile(unsigned samplingTimestepMultiple);
//
//    /**
//     * Find the solution of one of the PDEs at a point in space
//     *
//     * @param rPoint the position in space
//     * @param rVariable the dependent variable of the PDE whose solution you want to find
//     *
//     * @return the solution of the required PDE at the given point.
//     */
//    double GetPdeSolutionAtPoint(const c_vector<double,DIM>& rPoint, const std::string& rVariable);
//
//    /**
//     * @return the solution to the PDE at this time step.
//     *
//     * @param rName The name of the dependent variable for the PDE in the vector mPdeAndBcCollection.
//     * This defaults to an empty string in the case there is only one PDE.
//     */
//    virtual Vec GetPdeSolution(const std::string& rName = "");
//
//};
//
//#endif /*FIELDSOLVEREHANDLER_HPP_*/
