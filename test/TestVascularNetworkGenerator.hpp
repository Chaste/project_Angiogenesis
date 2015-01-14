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

#ifndef TESTVASCULARNETWORKGENERATOR_HPP_
#define TESTVASCULARNETWORKGENERATOR_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <algorithm>

#include "CellsGenerator.hpp"
#include "CaBasedCellPopulation.hpp"
#include "VolumeConstraintPottsUpdateRule.hpp"
#include "PottsMeshGenerator.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "DiffusionCaUpdateRule.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "ArchiveOpener.hpp"
#include "WildTypeCellMutationState.hpp"
#include "SmartPointers.hpp"
#include "CellLabel.hpp"
#include "FileComparison.hpp"

// Cell writers
#include "CellAgesWriter.hpp"
#include "CellAncestorWriter.hpp"
#include "CellProliferativePhasesWriter.hpp"
#include "CellVolumesWriter.hpp"
#include "CellMutationStatesWriter.hpp"

// Cell population writers
#include "CellMutationStatesCountWriter.hpp"
#include "CellProliferativeTypesCountWriter.hpp"
#include "CellProliferativePhasesCountWriter.hpp"

#include "VascularNetworkGenerator.hpp"

#include "PetscSetupAndFinalize.hpp"

class TestVascularNetworkGenerator : public AbstractCellBasedTestSuite
{
public:

    void TestVascularNetworkGeneratorWithHexagonallyTesselatedVesselNetwork() throw(Exception)
    {

        try
        {

            // Create a simple 2D PottsMesh
            unsigned width = 50;
            unsigned height = 25;
            PottsMeshGenerator<2> generator(width, 0, 0, height, 0, 0);
            PottsMesh<2>* p_mesh = generator.GetMesh();

            boost::shared_ptr<CaVessel<2> > prototypeVessel(new CaVessel<2>());
            prototypeVessel->GetCollectionOfIntraVascularChemicals().AddIntraVascularChemical("Oxygen",Concentration(0.0,"unitless"),3.0*pow(10.0,(-6)));
            VascularNetworkGenerator<2> vascularNetworkGenerator(prototypeVessel);

            boost::shared_ptr<CaVascularNetwork<2> > vascularNetwork = vascularNetworkGenerator.GenerateHexagonallyTesselatedVascularNetwork(p_mesh,width,height);

            // write vascular cell data out to file

            std::string output_directory = "TestVascularNetworkGenerator";
            std::string filename;
            std::stringstream sstm;
            OutputFileHandler output_file_handler(output_directory, false);

            sstm << output_file_handler.GetOutputDirectoryFullPath() << "HexagonallyTessellatedVesselNetwork.vtk";

            filename = sstm.str();

            vascularNetwork->SaveVasculatureDataToFile(filename);

        }
        catch(Exception &e)
        {
            std::cout << e.GetMessage() << std::endl;
        }

    }

};

#endif /*TESTVASCULARNETWORKGENERATOR_HPP_*/
