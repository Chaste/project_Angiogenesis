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
#include "OnLatticeSimulationWrapper.hpp"
#include "NodeBasedSimulationWrapper.hpp"
#include "AbstractCellPopulation.hpp"
#include "VascularTumourSolver.hpp"
#include "VascularTumourModifier.hpp"
#include "SimulationManager.hpp"
#include "LQRadiotherapyCellKiller.hpp"

using namespace boost::python;

// Make the module
BOOST_PYTHON_MODULE(_simulation)
{
    class_<SimulationManager, boost::shared_ptr<SimulationManager>, boost::noncopyable>("SimulationManager")
            .def("Setup", &SimulationManager::Setup)
            .def("TearDown", &SimulationManager::TearDown)
            .def("SetEndTimeAndNumberOfTimeSteps", &SimulationManager::SetEndTimeAndNumberOfTimeSteps)
    ;

    class_<OnLatticeSimulationWrapper, boost::shared_ptr<OnLatticeSimulationWrapper> >("OnLatticeSimulationWrapper")
            .def("SetCellPopulation", &OnLatticeSimulationWrapper::SetCellPopulation)
            .def("SetOutputDirectory", &OnLatticeSimulationWrapper::SetOutputDirectory)
            .def("SetDt", &OnLatticeSimulationWrapper::SetDt)
            .def("SetSamplingTimestepMultiple", &OnLatticeSimulationWrapper::SetSamplingTimestepMultiple)
            .def("SetEndTime", &OnLatticeSimulationWrapper::SetEndTime)
            .def("Solve", &OnLatticeSimulationWrapper::Solve)
            .def("GetOutputPopulations", &OnLatticeSimulationWrapper::GetOutputPopulations)
            .def("SetUseRadiotherapyCellKiller", &OnLatticeSimulationWrapper::SetUseRadiotherapyCellKiller)
            .def("SetNetwork", &OnLatticeSimulationWrapper::SetNetwork)
            .def("SetVesselDistanceTolerance", &OnLatticeSimulationWrapper::SetVesselDistanceTolerance)
            .def("SetRadiotherapyHitTimes", &OnLatticeSimulationWrapper::SetRadiotherapyHitTimes)
            .def("SetRadiotherapyDose", &OnLatticeSimulationWrapper::SetRadiotherapyDose)
            .def("SetOerAlphaMax", &OnLatticeSimulationWrapper::SetOerAlphaMax)
            .def("SetOerAlphaMin", &OnLatticeSimulationWrapper::SetOerAlphaMin)
            .def("SetOerBetaMax", &OnLatticeSimulationWrapper::SetOerBetaMax)
            .def("SetOerBetaMin", &OnLatticeSimulationWrapper::SetOerBetaMin)
            .def("SetOerConstant", &OnLatticeSimulationWrapper::SetOerConstant)
            .def("SetAlphaMax", &OnLatticeSimulationWrapper::SetAlphaMax)
            .def("SetBetaMax", &OnLatticeSimulationWrapper::SetBetaMax)
            .def("UseOer", &OnLatticeSimulationWrapper::UseOer)
    ;

    class_<VascularTumourModifier<3> >("VascularTumourModifier")
            .def("OutputSimulationModifierParameters", &VascularTumourModifier<3>::OutputSimulationModifierParameters)
            .def("UpdateCellData", &VascularTumourModifier<3>::UpdateCellData)
            .def("SetVascularTumourSolver", &VascularTumourModifier<3>::SetVascularTumourSolver)
            .def("SetCellDataUpdateLabels", &VascularTumourModifier<3>::SetCellDataUpdateLabels)
            .def("SetupSolve", &VascularTumourModifier<3>::SetupSolve)
            .def("UpdateAtEndOfTimeStep", &VascularTumourModifier<3>::UpdateAtEndOfTimeStep)
    ;

    class_<VascularTumourSolver<3> >("VascularTumourSolver")
            .def("SetRegressionSolver", &VascularTumourSolver<3>::SetRegressionSolver)
            .def("Run", &VascularTumourSolver<3>::Run)
            .def("Increment", &VascularTumourSolver<3>::Increment)
            .def("SetVesselNetwork", &VascularTumourSolver<3>::SetVesselNetwork)
            .def("UpdateCellData", &VascularTumourSolver<3>::UpdateCellData)
            .def("Setup", &VascularTumourSolver<3>::Setup)
            .def("SetupFromModifier", &VascularTumourSolver<3>::SetupFromModifier)
            .def("SetStructuralAdaptationSolver", &VascularTumourSolver<3>::SetStructuralAdaptationSolver)
            .def("SetOutputFrequency", &VascularTumourSolver<3>::SetOutputFrequency)
            .def("SetOutputFileHandler", &VascularTumourSolver<3>::SetOutputFileHandler)
            .def("SetAngiogenesisSolver", &VascularTumourSolver<3>::SetAngiogenesisSolver)
            .def("GetHybridSolvers", &VascularTumourSolver<3>::GetHybridSolvers)
            .def("AddHybridSolver", &VascularTumourSolver<3>::AddHybridSolver)
    ;
}
#endif // CHASTE_ANGIOGENESIS_PYTHON
