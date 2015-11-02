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

#include "AngiogenesisModifier.hpp"
#include "MeshBasedCellPopulation.hpp"
//#include "SimpleCellPopulation.hpp"
#include "boost/lexical_cast.hpp"
#include "Debug.hpp"

template<unsigned DIM>
AngiogenesisModifier<DIM>::AngiogenesisModifier()
    : AbstractCellBasedSimulationModifier<DIM>(),
      mpSolver(),
      mOutputDirectory()
{
}

template<unsigned DIM>
AngiogenesisModifier<DIM>::~AngiogenesisModifier()
{
}

template<unsigned DIM>
void AngiogenesisModifier<DIM>::SetAngiogenesisSolver(boost::shared_ptr<AbstractAngiogenesisSolver<DIM> > pSolver)
{
    mpSolver = pSolver;
}

template<unsigned DIM>
void AngiogenesisModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // If there is an angiogenesis solver solve for the upcoming step
    if(mpSolver)
    {
        // Increment the solver
        mpSolver->Increment();
    }

    // Update the cell data
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void AngiogenesisModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    mOutputDirectory = outputDirectory;

    // If there is an angiogenesis solver solve for the upcoming step
    if(mpSolver)
    {
        for(unsigned idx=0; idx<mpSolver->GetPdeSolvers().size(); idx++)
        {
            for(unsigned jdx=0; jdx<mpSolver->GetPdeSolvers()[idx]->GetPde()->GetDiscreteSources().size(); jdx++)
            {
                if(mpSolver->GetPdeSolvers()[idx]->GetPde()->GetDiscreteSources()[jdx]->GetType() == SourceType::CELL)
                {
                    mpSolver->GetPdeSolvers()[idx]->GetPde()->GetDiscreteSources()[jdx]->SetCellPopulation(rCellPopulation);
                }
            }
        }

        // Set the output directory
        OutputFileHandler output_file_handler(mOutputDirectory, false);
        mpSolver->SetOutputDirectory(output_file_handler.GetOutputDirectoryFullPath());
        mpSolver->Increment();
    }

    // Update the cell data
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void AngiogenesisModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Make sure the cell population is updated
    rCellPopulation.Update();

    if(mpSolver)
    {
        // Get the cell locations
        std::vector<c_vector<double, DIM> > locations;
        std::vector<CellPtr> cell_vector;
        for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
             cell_iter != rCellPopulation.End(); ++cell_iter)
        {
            locations.push_back(rCellPopulation.GetLocationOfCellCentre(*cell_iter));
            cell_vector.push_back(*cell_iter);
        }

        // Get the pde solutions at each cell location
        for(unsigned idx=0; idx<mpSolver->GetPdeSolvers().size(); idx++)
        {
            std::string label = mpSolver->GetPdeSolvers()[idx]->GetPde()->GetVariableName();
            std::vector<double> sampled_solution = mpSolver->GetPdeSolvers()[idx]->GetSolutionAtPoints(locations, label);

            if(label == "oxygen")
            {
                for(unsigned jdx=0;jdx<cell_vector.size();jdx++)
                {
                    cell_vector[jdx]->GetCellData()->SetItem("oxygen", sampled_solution[jdx]);
                }
            }
            else
            {
                for(unsigned jdx=0;jdx<cell_vector.size();jdx++)
                {
                    cell_vector[jdx]->GetCellData()->SetItem(label, sampled_solution[jdx]);
                }
            }
        }
    }
    else
    {
        for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
             cell_iter != rCellPopulation.End();
             ++cell_iter)
        {
            cell_iter->GetCellData()->SetItem("empty", 0.0);
        }
    }
}

template<unsigned DIM>
void AngiogenesisModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
//template class AngiogenesisModifier<1>;
template class AngiogenesisModifier<2>;
template class AngiogenesisModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS1(AngiogenesisModifier, 2)
EXPORT_TEMPLATE_CLASS1(AngiogenesisModifier, 3)
