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
#include <vector>
#include <boost/python.hpp>
#include <boost/python/module.hpp>
#include <boost/python/stl_iterator.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#define _BACKWARD_BACKWARD_WARNING_H 1 //Cut out the vtk deprecated warning for now (gcc4.3)
#include <vtkSmartPointer.h>
#include <vtkImageData.h>
#include <vtkPolyData.h>
#include "FiniteDifferenceSolver.hpp"
#include "FiniteElementSolver.hpp"
#include "DistanceMap.hpp"
#include "FunctionMap.hpp"
#include "AbstractRegularGridHybridSolver.hpp"
#include "HybridBoundaryCondition.hpp"
#include "GreensFunctionSolver.hpp"
#include "CellStateDependentDiscreteSource.hpp"
#include "PdeWrappers.hpp"

#include "converters.hpp"

using namespace boost::python;


// Make the module
BOOST_PYTHON_MODULE(_pde)
{
    class_<AbstractHybridSolverWrap, boost::shared_ptr<AbstractHybridSolverWrap>, boost::noncopyable>("AbstractHybridSolver")
        .def("Solve", pure_virtual(&AbstractHybridSolver<3>::Solve))
        .def("Update", pure_virtual(&AbstractHybridSolver<3>::Update))
        .def("Write", pure_virtual(&AbstractHybridSolver<3>::Write))
        .def("Setup", pure_virtual(&AbstractHybridSolver<3>::Setup))
        .def("UpdateCellData", pure_virtual(&AbstractHybridSolver<3>::UpdateCellData))
        .def("SetLabel", &AbstractHybridSolver<3>::SetLabel)
        .def("SetFileName", &AbstractHybridSolver<3>::SetFileName)
        .def("SetFileHandler", &AbstractHybridSolver<3>::SetFileHandler)
        .def("SetVesselNetwork", &AbstractHybridSolver<3>::SetVesselNetwork)
        .def("SetCellPopulation", &AbstractHybridSolver<3>::SetCellPopulation)
        .def("AddBoundaryCondition", &AbstractHybridSolver<3>::AddBoundaryCondition)
        .def("GetPde", &AbstractHybridSolver<3>::GetPde)
        .def("GetNonLinearPde", &AbstractHybridSolver<3>::GetNonLinearPde)
        .def("SetPde", &AbstractHybridSolver<3>::SetPde)
        .def("SetNonLinearPde", &AbstractHybridSolver<3>::SetNonLinearPde)
        .def("SetWriteSolution", &AbstractHybridSolver<3>::SetWriteSolution)
        .def("GetSolutionAtPoints", pure_virtual(&AbstractHybridSolver<3>::GetSolutionAtPoints))
        .def("GetSolutionAtGridPoints", pure_virtual(&AbstractHybridSolver<3>::GetSolutionAtGridPoints))
        ;

    class_<AbstractRegularGridHybridSolverWrap, boost::shared_ptr<AbstractRegularGridHybridSolverWrap>, boost::noncopyable, bases<AbstractHybridSolver<3> >  >("AbstractRegularGridHybridSolver")
        .def("Solve", pure_virtual(&AbstractRegularGridHybridSolver<3>::Solve))
        .def("Update", pure_virtual(&AbstractRegularGridHybridSolver<3>::Update))
        .def("SetGrid", &AbstractRegularGridHybridSolver<3>::SetGrid)
        .def("GetGrid", &AbstractRegularGridHybridSolver<3>::GetGrid)
        .def("GetVtkSolution", &AbstractRegularGridHybridSolver<3>::GetVtkSolution, &AbstractRegularGridHybridSolverWrap::default_GetVtkSolution)
        .def("GetPointSolution", &AbstractRegularGridHybridSolver<3>::GetPointSolution)
        .def("GetSolutionAtPoints", &AbstractRegularGridHybridSolver<3>::GetSolutionAtPoints)
        .def("GetSolutionAtGridPoints", &AbstractRegularGridHybridSolver<3>::GetSolutionAtGridPoints)
        ;

    class_<DistanceMap<3>, boost::shared_ptr<DistanceMap<3u> >, bases<AbstractRegularGridHybridSolver<3> > >("DistanceMap")
        .def("Solve", &DistanceMap<3>::Solve)
    ;

    class_<FunctionMap<3>, boost::shared_ptr<FunctionMap<3u> >, bases<AbstractRegularGridHybridSolver<3> >  >("FunctionMap")
        .def("Solve", &FunctionMap<3>::Solve)
        .def("UpdateSolution", &FunctionMap<3>::UpdateSolution)
    ;

    class_<FiniteDifferenceSolver<3>, boost::shared_ptr<FiniteDifferenceSolver<3u> >, bases<AbstractRegularGridHybridSolver<3> > >("FiniteDifferenceSolver")
        .def("Solve", &FiniteDifferenceSolver<3>::Solve)
        .def("Update", &FiniteDifferenceSolver<3>::Update)
        .def("UpdateBoundaryConditionsEachSolve", &FiniteDifferenceSolver<3>::UpdateBoundaryConditionsEachSolve)
    ;

    class_<GreensFunctionSolver<3>, boost::shared_ptr<GreensFunctionSolver<3u> >, bases<AbstractRegularGridHybridSolver<3> > >("GreensFunctionSolver")
        .def("Solve", &GreensFunctionSolver<3>::Solve)
        .def("SetSubSegmentCutoff", &GreensFunctionSolver<3>::SetSubSegmentCutoff)
        .def("Setup", &AbstractRegularGridHybridSolver<3>::Setup)
    ;

    class_<FiniteElementSolver<3>, boost::shared_ptr<FiniteElementSolver<3u> > >("FiniteElementSolver")
        .def("Solve", &FiniteElementSolver<3>::Solve)
        .def("Update", &FiniteElementSolver<3>::Update)
        .def("SetMesh", &FiniteElementSolver<3>::SetMesh)
        .def("GetSolutionAtPoints", &FiniteElementSolver<3>::GetSolutionAtPoints)
        .def("GetSolutionAtGridPoints", &AbstractRegularGridHybridSolver<3>::GetSolutionAtGridPoints)
    ;

    enum_<BoundaryConditionType::Value>("BoundaryConditionType")
        .value("FACET", BoundaryConditionType::FACET)
        .value("OUTER", BoundaryConditionType::OUTER)
        .value("VESSEL_LINE", BoundaryConditionType::VESSEL_LINE)
        .value("VESSEL_VOLUME", BoundaryConditionType::VESSEL_VOLUME)
        .value("CELL", BoundaryConditionType::CELL)
        .value("IN_PART", BoundaryConditionType::IN_PART)
        .value("POINT", BoundaryConditionType::POINT)
        ;

    enum_<BoundaryConditionSource::Value>("BoundaryConditionSource")
        .value("LABEL_BASED", BoundaryConditionSource::LABEL_BASED)
        .value("PRESCRIBED", BoundaryConditionSource::PRESCRIBED)
        ;

    class_<HybridBoundaryCondition<3>, boost::shared_ptr<HybridBoundaryCondition<3u> > >("HybridBoundaryCondition")
        .def("SetValue", &HybridBoundaryCondition<3>::SetValue)
        .def("SetType", &HybridBoundaryCondition<3>::SetType)
        .def("SetSource", &HybridBoundaryCondition<3>::SetSource)
        .def("SetLabelName", &HybridBoundaryCondition<3>::SetLabelName)
        .def("SetDomain", &HybridBoundaryCondition<3>::SetDomain)
        .def("SetNetwork", &HybridBoundaryCondition<3>::SetNetwork)
        .def("SetPoints", &HybridBoundaryCondition<3>::SetPoints)
    ;

    enum_<SourceType::Value>("DiscreteSourceType")
        .value("POINT", SourceType::POINT)
        .value("VESSEL", SourceType::VESSEL)
        .value("CELL", SourceType::CELL)
        .value("SOLUTION", SourceType::SOLUTION)
        ;

    enum_<SourceStrength::Value>("DiscreteSourceStrength")
        .value("LABEL", SourceStrength::LABEL)
        .value("PRESCRIBED", SourceStrength::PRESCRIBED)
        ;

    class_<DiscreteSource<3>, boost::shared_ptr<DiscreteSource<3u> > >("DiscreteSource")
        .def("SetSolution", &DiscreteSource<3>::SetSolution)
        .def("SetType", &DiscreteSource<3>::SetType)
        .def("SetSource", &DiscreteSource<3>::SetSource)
        .def("SetLabelName", &DiscreteSource<3>::SetLabelName)
        .def("SetRegularGrid", &DiscreteSource<3>::SetRegularGrid)
        .def("SetMesh", &DiscreteSource<3>::SetMesh)
        .def("SetValue", &DiscreteSource<3>::SetValue)
        .def("SetIsLinearInSolution", &DiscreteSource<3>::SetIsLinearInSolution)
    ;

    class_<CellStateDependentDiscreteSource<3>, boost::shared_ptr<CellStateDependentDiscreteSource<3u> >, bases<DiscreteSource<3> > >("CellStateDependentDiscreteSource")
        .def("SetSolution", &CellStateDependentDiscreteSource<3>::SetSolution)
        .def("SetType", &CellStateDependentDiscreteSource<3>::SetType)
        .def("SetSource", &CellStateDependentDiscreteSource<3>::SetSource)
        .def("SetLabelName", &CellStateDependentDiscreteSource<3>::SetLabelName)
        .def("SetRegularGrid", &CellStateDependentDiscreteSource<3>::SetRegularGrid)
        .def("SetMesh", &CellStateDependentDiscreteSource<3>::SetMesh)
        .def("SetValue", &CellStateDependentDiscreteSource<3>::SetValue)
        .def("SetIsLinearInSolution", &CellStateDependentDiscreteSource<3>::SetIsLinearInSolution)
        .def("SetStateRateMap", &CellStateDependentDiscreteSource<3>::SetStateRateMap)
    ;

    class_<HybridLinearEllipticPde<3>, boost::shared_ptr<HybridLinearEllipticPde<3u> > >("HybridLinearEllipticPde")
        .def("SetIsotropicDiffusionConstant", &HybridLinearEllipticPde<3>::SetIsotropicDiffusionConstant)
        .def("SetContinuumLinearInUTerm", &HybridLinearEllipticPde<3>::SetContinuumLinearInUTerm)
        .def("SetContinuumConstantInUTerm", &HybridLinearEllipticPde<3>::SetContinuumConstantInUTerm)
        .def("AddDiscreteSource", &HybridLinearEllipticPde<3>::AddDiscreteSource)
        .def("UpdateDiscreteSourceStrengths", &HybridLinearEllipticPde<3>::UpdateDiscreteSourceStrengths)
        .def("SetRegularGrid", &HybridLinearEllipticPde<3>::SetRegularGrid)
        .def("SetVariableName", &HybridLinearEllipticPde<3>::SetVariableName)
    ;

    class_<HybridNonLinearEllipticPde<3>, boost::shared_ptr<HybridNonLinearEllipticPde<3u> > >("HybridNonLinearEllipticPde")
        .def("SetIsotropicDiffusionConstant", &HybridNonLinearEllipticPde<3>::SetIsotropicDiffusionConstant)
        .def("SetContinuumLinearInUTerm", &HybridNonLinearEllipticPde<3>::SetContinuumLinearInUTerm)
        .def("SetContinuumConstantInUTerm", &HybridNonLinearEllipticPde<3>::SetContinuumConstantInUTerm)
        .def("AddDiscreteSource", &HybridNonLinearEllipticPde<3>::AddDiscreteSource)
        .def("UpdateDiscreteSourceStrengths", &HybridNonLinearEllipticPde<3>::UpdateDiscreteSourceStrengths)
        .def("SetRegularGrid", &HybridNonLinearEllipticPde<3>::SetRegularGrid)
        .def("SetThreshold", &HybridNonLinearEllipticPde<3>::SetThreshold)
    ;
}
#endif // CHASTE_ANGIOGENESIS_PYTHON
