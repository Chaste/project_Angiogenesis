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

#include <algorithm>
#include <vector>
#include <iostream>
#include "StructuralAdaptationSolver.hpp"
#include "SimulationTime.hpp"

template<unsigned DIM>
StructuralAdaptationSolver<DIM>::StructuralAdaptationSolver()
	: mTolerance(0.0001),
	  mTimeIncrement(0.0001),
	  mWriteOutput(false),
	  mOutputFileName("StructuralAdaptationSolverProgress"),
	  mMaxIterations(100000)
{

}

template<unsigned DIM>
StructuralAdaptationSolver<DIM>::~StructuralAdaptationSolver()
{

}

template<unsigned DIM>
double StructuralAdaptationSolver<DIM>::GetTolerance() const
{
    return mTolerance;
}

template<unsigned DIM>
bool StructuralAdaptationSolver<DIM>::GetWriteOutput() const
{
    return mWriteOutput;
}


template<unsigned DIM>
std::string StructuralAdaptationSolver<DIM>::GetOutputFileName() const
{
    return mOutputFileName;
}

template<unsigned DIM>
double StructuralAdaptationSolver<DIM>::GetTimeIncrement() const
{
    return mTimeIncrement;
}


template<unsigned DIM>
void StructuralAdaptationSolver<DIM>::SetTolerance(double tolerance)
{
	mTolerance = tolerance;
}


template<unsigned DIM>
void StructuralAdaptationSolver<DIM>::SetTimeIncrement(double timeIncrement)
{
	mTimeIncrement = timeIncrement;
}

template<unsigned DIM>
void StructuralAdaptationSolver<DIM>::SetMaxIterations(unsigned iterations)
{
	mMaxIterations = iterations;
}

template<unsigned DIM>
void StructuralAdaptationSolver<DIM>::SetWriteOutput(bool writeFlag)
{
	mWriteOutput = writeFlag;
}


template<unsigned DIM>
void StructuralAdaptationSolver<DIM>::SetOutputFileName(std::string filename)
{
	mOutputFileName = filename;
}


template<unsigned DIM>
void StructuralAdaptationSolver<DIM>::Implement(boost::shared_ptr<CaVascularNetwork<DIM> > vascularNetwork)
{

    std::ofstream out;

    double max_radius_relative_change = 1.0;
    unsigned iteration = 0;
    double time = 0.0;

    if (mWriteOutput)
    {
        out.open(mOutputFileName.c_str());
    }

    if (out.is_open())
    {
        out << "\n";
        out << "#Iteration   Maximum relative change in radius in network\n\n";
    }

    std::vector<double> previous_radii;
    std::vector<boost::shared_ptr<CaVesselSegment<DIM> > > segments = vascularNetwork->GetVesselSegments();

    for (unsigned segment_index = 0; segment_index < segments.size(); segment_index++)
    {
    	previous_radii.push_back(segments[segment_index]->GetRadius());
    }

    while (max_radius_relative_change > mTolerance && time < SimulationTime::Instance()->GetTimeStep()*60.0 && iteration < mMaxIterations)
    {
        time += mTimeIncrement;
        iteration++;

        Iterate(vascularNetwork);

        std::vector<double> relative_change;

        for (unsigned segment_index = 0; segment_index < segments.size(); segment_index++)
        {
        	double current_radius = segments[segment_index]->GetRadius();

        	relative_change.push_back(fabs(1.0 - current_radius/previous_radii[segment_index]));
        	previous_radii[segment_index] = current_radius;
        }

        max_radius_relative_change = *( std::max_element(relative_change.begin(), relative_change.end()));

        if (out.is_open())
        {
            out<< std::setw(6) << iteration << std::setw(20) << max_radius_relative_change<< "\n";
        }

    }

    if (out.is_open())
    {
        out.close();
    }

}

template<unsigned DIM>
void StructuralAdaptationSolver<DIM>::WriteToFile(std::string parameterFileName)
{
    std::ofstream out;

    out.open(parameterFileName.c_str(), std::ios::app); // append to file
    out << "\nModule: StructuralAdaptationSolver\n";
    out << "\n---------------\n";
    out.close();
}

// Explicit instantiation
template class StructuralAdaptationSolver<2>;
template class StructuralAdaptationSolver<3>;
