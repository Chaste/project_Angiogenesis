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

#include <vector>
#include <boost/python.hpp>
#include <boost/python/module.hpp>
#include "AngiogenesisSolver.hpp"
#include "Owen2011MigrationRule.hpp"
#include "Owen2011SproutingRule.hpp"
#include "OffLatticeMigrationRule.hpp"
#include "OffLatticeSproutingRule.hpp"
#include "AbstractMigrationRule.hpp"
#include "AbstractSproutingRule.hpp"

using namespace boost::python;

// Make the module
BOOST_PYTHON_MODULE(_angiogenesis)
{

    class_<AngiogenesisSolver<3> >("AngiogenesisSolver")
            .def("SetVesselNetwork", &AngiogenesisSolver<3>::SetVesselNetwork)
            .def("SetVesselGrid", &AngiogenesisSolver<3>::SetVesselGrid)
            .def("SetMigrationRule", &AngiogenesisSolver<3>::SetMigrationRule)
            .def("SetSproutingRule", &AngiogenesisSolver<3>::SetSproutingRule)
            .def("SetOutputFileHandler", &AngiogenesisSolver<3>::SetOutputFileHandler)
            .def("Run", &AngiogenesisSolver<3>::Run)
    ;

    class_<AbstractMigrationRule<3>, boost::shared_ptr<AbstractMigrationRule<3u> > >("AbstractMigrationRule")
            .def("SetHybridSolver", &AbstractMigrationRule<3>::SetHybridSolver)
            .def("SetNetwork", &AbstractMigrationRule<3>::SetNetwork)
    ;

    class_<Owen2011MigrationRule<3>, boost::shared_ptr<Owen2011MigrationRule<3u> >, bases<AbstractMigrationRule<3> > >("Owen2011MigrationRule")
            .def("SetGrid", &Owen2011MigrationRule<3>::SetGrid)
            .def("SetHybridSolver", &Owen2011MigrationRule<3>::SetHybridSolver)
            .def("SetCellMotilityParameter", &Owen2011MigrationRule<3>::SetCellMotilityParameter)
            .def("SetCellChemotacticParameter", &Owen2011MigrationRule<3>::SetCellChemotacticParameter)
            .def("SetNetwork", &Owen2011MigrationRule<3>::SetNetwork)
    ;

    class_<AbstractSproutingRule<3>, boost::shared_ptr<AbstractSproutingRule<3u> > >("AbstractSproutingRule")
            .def("SetHybridSolver", &AbstractSproutingRule<3>::SetHybridSolver)
            .def("SetSproutingProbability", &AbstractSproutingRule<3>::SetSproutingProbability)
            .def("SetGrid", &AbstractSproutingRule<3>::SetGrid)
            .def("SetVesselNetwork", &AbstractSproutingRule<3>::SetVesselNetwork)
    ;

    class_<Owen2011SproutingRule<3>, boost::shared_ptr<Owen2011SproutingRule<3u> >, bases<AbstractSproutingRule<3> > >("Owen2011SproutingRule")
            .def("SetHybridSolver", &Owen2011SproutingRule<3>::SetHybridSolver)
            .def("SetSproutingProbability", &Owen2011SproutingRule<3>::SetSproutingProbability)
            .def("SetGrid", &Owen2011SproutingRule<3>::SetGrid)
            .def("SetVesselNetwork", &Owen2011SproutingRule<3>::SetVesselNetwork)
    ;

    class_<OffLatticeSproutingRule<3>, boost::shared_ptr<OffLatticeSproutingRule<3u> >, bases<AbstractSproutingRule<3> > >("OffLatticeSproutingRule")
            .def("SetHybridSolver", &OffLatticeSproutingRule<3>::SetHybridSolver)
            .def("SetSproutingProbability", &OffLatticeSproutingRule<3>::SetSproutingProbability)
            .def("SetVesselNetwork", &OffLatticeSproutingRule<3>::SetVesselNetwork)
    ;

    class_<OffLatticeMigrationRule<3>, boost::shared_ptr<OffLatticeMigrationRule<3u> >, bases<AbstractMigrationRule<3> > >("OffLatticeMigrationRule")
            .def("SetHybridSolver", &OffLatticeMigrationRule<3>::SetHybridSolver)
            .def("SetNetwork", &OffLatticeMigrationRule<3>::SetNetwork)
            .def("SetSproutingVelocity", &OffLatticeMigrationRule<3>::SetSproutingVelocity)
            .def("SetChemotacticStrength", &OffLatticeMigrationRule<3>::SetChemotacticStrength)
            .def("SetAttractionStrength", &OffLatticeMigrationRule<3>::SetAttractionStrength)
    ;
}
