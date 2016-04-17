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

#include <numeric>
#include <cstdlib>
#define _BACKWARD_BACKWARD_WARNING_H 1 //Cut out the vtk deprecated warning for now (gcc4.3)
#include <vtkXMLImageDataWriter.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkImageData.h>
#include <vtkSmartPointer.h>
#include <boost/lexical_cast.hpp>
#include "CaVessel.hpp"
#include "CaVesselSegment.hpp"
#include "ChastePoint.hpp"
#include "LinearSystem.hpp"
#include "ReplicatableVector.hpp"
#include "UblasMatrixInclude.hpp"

#include "GreensFunctionSolver.hpp"

template<unsigned DIM>
GreensFunctionSolver<DIM>::GreensFunctionSolver()
    : AbstractRegularGridHybridSolver<DIM>(),
      mpDomain(),
      mSinkCoordinates(),
      mSinkPointMap(),
      mSubSegmentCoordinates(),
      mSubSegmentLengths(),
      mSinkRates(),
      mSourceRates(),
      mSegmentConcentration(),
      mTissueConcentration(),
      mSegmentPointMap(),
      mGtt(),
      mGvv(),
      mGvt(),
      mGtv(),
      mSubsegmentCutoff(1.0)
{

}

template<unsigned DIM>
GreensFunctionSolver<DIM>::~GreensFunctionSolver()
{

}

template<unsigned DIM>
void GreensFunctionSolver<DIM>::Solve()
{
    // Set up the sub-segment and tissue point co-ordinates
    GenerateSubSegments();
    GenerateTissuePoints();

    // Generate the greens function matrices
    UpdateGreensFunctionMatrices(1, 1, 1, 1);

    // Get the sink rates
    unsigned number_of_sinks = mSinkCoordinates.size();
    double sink_rate = this->mpPde->ComputeConstantInUSourceTerm();
    double sink_volume = pow(this->mpRegularGrid->GetSpacing(), 3);
    mSinkRates = std::vector<double>(number_of_sinks, sink_rate * sink_volume);
    double total_sink_rate = std::accumulate(mSinkRates.begin(), mSinkRates.end(), 0.0);

    // Get the sink substance demand on each vessel subsegment
    unsigned number_of_subsegments = mSubSegmentCoordinates.size();
    double diffusivity = this->mpPde->ComputeIsotropicDiffusionTerm();
    std::vector<double> sink_demand_per_subsegment(number_of_subsegments, 0.0);
    for (unsigned idx = 0; idx < number_of_subsegments; idx++)
    {
        for (unsigned jdx = 0; jdx < number_of_sinks; jdx++)
        {
            sink_demand_per_subsegment[idx] += ((*mGvt)[idx][jdx] / diffusivity) * mSinkRates[jdx];
        }
    }

    mSegmentConcentration = std::vector<double>(number_of_subsegments, 1.0);
    mTissueConcentration = std::vector<double>(number_of_sinks, 0.0);

    // Solve for the subsegment source rates required to meet the sink substance demand
    double tolerance = 1.e-10;
    double g0 = 0.0;
    mSourceRates = std::vector<double>(number_of_subsegments, 0.0);

    LinearSystem linear_system(number_of_subsegments + 1, number_of_subsegments + 1);
    linear_system.SetKspType("bcgs");

    for (unsigned iteration = 0; iteration < 10; iteration++)
    {
        linear_system.AssembleIntermediateLinearSystem();
        for (unsigned i = 0; i < number_of_subsegments; i++)
        {
            linear_system.SetRhsVectorElement(i, mSegmentConcentration[i] - sink_demand_per_subsegment[i]);
        }

        linear_system.SetRhsVectorElement(number_of_subsegments, -total_sink_rate);

        // Set up Linear system matrix
        for (unsigned iter = 0; iter < number_of_subsegments; iter++)
        {
            for (unsigned jter = 0; jter < number_of_subsegments; jter++)
            {
                linear_system.SetMatrixElement(iter, jter, (*mGvv)[iter][jter] / diffusivity);
            }
            linear_system.SetMatrixElement(number_of_subsegments, iter, 1.0);
            linear_system.SetMatrixElement(iter, number_of_subsegments, 1.0);
        }
        linear_system.SetMatrixElement(number_of_subsegments, number_of_subsegments, 0.0);

        // Solve the linear system
        linear_system.AssembleFinalLinearSystem();
        ReplicatableVector soln_repl(linear_system.Solve());

        // Populate the solution vector
        std::vector<double> solution_vector(number_of_subsegments + 1);
        for (unsigned row = 0; row < number_of_subsegments + 1; row++)
        {
            (solution_vector)[row] = soln_repl[row];
        }

        // Check convergence
        bool all_in_tolerance = true;
        for (unsigned i = 0; i < number_of_subsegments; i++)
        {
            double diff = std::abs(mSourceRates[i] - solution_vector[i]);
            if (diff > tolerance)
            {
                all_in_tolerance = false;
                break;
            }
        }
        // Retrieve the solution
        for (unsigned i = 0; i < number_of_subsegments; i++)
        {
            mSourceRates[i] = solution_vector[i];
        }
        g0 = solution_vector[number_of_subsegments];

        if (all_in_tolerance)
        {
            break;
        }
        else
        {
            if (iteration == 9)
            {
                std::cout << "Did not converge\n";
            }
        }
    }

    // Get the tissue concentration and write the solution
    for (unsigned i = 0; i < number_of_sinks; i++)
    {
        mTissueConcentration[i] = 0.0;
        for (unsigned j = 0; j < number_of_sinks; j++)
        {
            mTissueConcentration[i] += (*mGtt)[i][j] * mSinkRates[j] / diffusivity;
        }

        for (unsigned j = 0; j < number_of_subsegments; j++)
        {
            mTissueConcentration[i] += (*mGtv)[i][j] * mSourceRates[j] / diffusivity;
        }
        mTissueConcentration[i] += g0;
    }

    std::map<std::string, std::vector<double> > segmentPointData;
    std::map<std::string, std::vector<double> > tissuePointData;
    tissuePointData[this->mpPde->GetVariableName()] = mTissueConcentration;

    this->mPointSolution = std::vector<double>(mTissueConcentration.size(), 0.0);
    for(unsigned idx=0; idx<mTissueConcentration.size(); idx++)
    {
        this->mPointSolution[idx] = mTissueConcentration[idx];
    }


    segmentPointData[this->mpPde->GetVariableName()] = mSegmentConcentration;
    tissuePointData["Sink Rate"] = mSinkRates;
    segmentPointData["Source Rate"] = mSourceRates;
    this->UpdateSolution(tissuePointData);
    if(this->mWriteSolution)
    {
        this->WriteSolution(segmentPointData);
    }
}

template<unsigned DIM>
void GreensFunctionSolver<DIM>::GenerateSubSegments()
{
    // Set up the sub-segment points and map to original segments
    double max_subsegment_length = mSubsegmentCutoff;

    std::vector<boost::shared_ptr<CaVessel<DIM> > > vessels = this->mpNetwork->GetVessels();
    typename std::vector<boost::shared_ptr<CaVessel<DIM> > >::iterator vessel_iter;
    typename std::vector<boost::shared_ptr<CaVesselSegment<DIM> > >::iterator segment_iter;

    // Iterate over all segments and store midpoints and lengths of subsegment regions for
    // the greens functions calculation. Create a map of subsegment index to the parent segment
    // for later use.
    for (vessel_iter = vessels.begin(); vessel_iter != vessels.end(); vessel_iter++)
    {
        std::vector<boost::shared_ptr<CaVesselSegment<DIM> > > segments = (*vessel_iter)->GetSegments();
        for (segment_iter = segments.begin(); segment_iter != segments.end(); segment_iter++)
        {
            double segment_length = (*segment_iter)->GetLength();

            // If the segment is shorter than the max length just use its mid-point
            if (segment_length < 1.01 * max_subsegment_length)
            {
                mSubSegmentCoordinates.push_back((*segment_iter)->GetMidPoint());
                mSubSegmentLengths.push_back(segment_length);
                mSegmentPointMap[mSubSegmentCoordinates.size() - 1] = (*segment_iter);
            }
            // Otherwise generate subsegment points along its length
            else
            {
                c_vector<double, DIM> start_point = (*segment_iter)->GetNode(0)->GetLocationVector();
                c_vector<double, DIM> end_point = (*segment_iter)->GetNode(1)->GetLocationVector();

                double subsegment_length = segment_length / max_subsegment_length;
                unsigned num_subsegments = std::floor(subsegment_length) + 1;
                subsegment_length = segment_length / double(num_subsegments);
                c_vector<double, DIM> increment = (end_point - start_point) / segment_length;

                for (unsigned i = 0; i < num_subsegments; i++)
                {
                    c_vector<double, DIM> location = start_point + (double(i) + 0.5) * increment * subsegment_length;
                    mSubSegmentCoordinates.push_back(ChastePoint<DIM>(location));
                    mSubSegmentLengths.push_back(subsegment_length);
                    mSegmentPointMap[mSubSegmentCoordinates.size() - 1] = (*segment_iter);
                }
            }
        }
    }
}

template<unsigned DIM>
void GreensFunctionSolver<DIM>::GenerateTissuePoints()
{
    unsigned num_points = this->mpRegularGrid->GetNumberOfPoints();
    mSinkCoordinates = std::vector<ChastePoint<DIM> >(num_points);
    mSinkPointMap = std::vector<unsigned>(num_points);
    for(unsigned idx=0; idx<num_points; idx++)
    {
        mSinkCoordinates[idx] = this->mpRegularGrid->GetLocationOf1dIndex(idx);
        mSinkPointMap[idx] = idx;
    }
}

template<unsigned DIM>
void GreensFunctionSolver<DIM>::UpdateGreensFunctionMatrices(bool updateGtt, bool updateGvv, bool updateGtv,
                                                                 bool updateGvt)
{
    // Get the Greens Function coefficient matrices
    if (updateGtt)
    {
        mGtt = GetTissueTissueInteractionMatrix();
    }
    if (updateGvv)
    {
        mGvv = GetVesselVesselInteractionMatrix();
    }
    if (updateGtv)
    {
        mGtv = GetTissueVesselInteractionMatrix();
    }
    if (updateGvt)
    {
        mGvt = GetVesselTissueInteractionMatrix();
    }
}

template<unsigned DIM>
boost::shared_ptr<boost::multi_array<double, 2> > GreensFunctionSolver<DIM>::GetVesselVesselInteractionMatrix()
{
    typedef boost::multi_array<double, 2>::index index;
    unsigned num_sub_segments = mSubSegmentCoordinates.size();
    double coefficient = 1.0 / (4.0 * M_PI);

    boost::shared_ptr<boost::multi_array<double, 2> > p_interaction_matrix(
            new boost::multi_array<double, 2>(boost::extents[num_sub_segments][num_sub_segments]));
    for (index iter = 0; iter < num_sub_segments; iter++)
    {
        for (index iter2 = 0; iter2 < num_sub_segments; iter2++)
        {
            if (iter <= iter2)
            {
                double distance = norm_2(mSubSegmentCoordinates[iter2].rGetLocation() - mSubSegmentCoordinates[iter].rGetLocation());
                double term;
                if (distance < mSegmentPointMap[iter]->GetRadius())
                {
                    double radius = mSegmentPointMap[iter]->GetRadius();
                    double max_segment_length = std::max(mSubSegmentLengths[iter], mSubSegmentLengths[iter2]);
                    double green_correction = 0.6 * std::exp(-0.45 * max_segment_length /radius);

                    if (iter != iter2)
                    {
                        distance = radius;
                    }

                    term = (1.298 / (1.0 + 0.297 * pow(max_segment_length /radius, 0.838))- green_correction * pow(distance / radius, 2)) *
                            coefficient/ radius;
                }
                else
                {
                    term = coefficient / distance;
                }
                (*p_interaction_matrix)[iter][iter2] = term;
                (*p_interaction_matrix)[iter2][iter] = term;
            }
        }
    }
    return p_interaction_matrix;
}

template<unsigned DIM>
boost::shared_ptr<boost::multi_array<double, 2> > GreensFunctionSolver<DIM>::GetTissueTissueInteractionMatrix()
{
    typedef boost::multi_array<double, 2>::index index;
    unsigned num_points = mSinkCoordinates.size();
    double coefficient = 1.0 / (4.0 * M_PI);
    double tissue_point_volume = pow(this->mpRegularGrid->GetSpacing(), 3);
    double equivalent_tissue_point_radius = pow(tissue_point_volume * 0.75 / M_PI, 0.333333);

    boost::shared_ptr<boost::multi_array<double, 2> > p_interaction_matrix(new boost::multi_array<double, 2>(boost::extents[num_points][num_points]));
    for (index iter = 0; iter < num_points; iter++)
    {
        for (index iter2 = 0; iter2 < num_points; iter2++)
        {
            if (iter < iter2)
            {
                double distance = norm_2(mSinkCoordinates[iter2].rGetLocation() - mSinkCoordinates[iter].rGetLocation());
                (*p_interaction_matrix)[iter][iter2] = coefficient / distance;
                (*p_interaction_matrix)[iter2][iter] = coefficient / distance;
            }
            else if (iter == iter2)
            {
                (*p_interaction_matrix)[iter][iter2] = 1.2 * coefficient / equivalent_tissue_point_radius;
            }
        }
    }
    return p_interaction_matrix;
}

template<unsigned DIM>
boost::shared_ptr<boost::multi_array<double, 2> > GreensFunctionSolver<DIM>::GetTissueVesselInteractionMatrix()
{
    typedef boost::multi_array<double, 2>::index index;
    unsigned num_sinks = mSinkCoordinates.size();
    unsigned num_subsegments = mSubSegmentCoordinates.size();

    double tissue_point_volume = pow(this->mpRegularGrid->GetSpacing(), 3);
    double equivalent_tissue_point_radius = pow(tissue_point_volume * 0.75 / M_PI, 0.333333);
    double coefficient = 1.0 / (4.0 * M_PI);

    boost::shared_ptr<boost::multi_array<double, 2> > p_interaction_matrix(new boost::multi_array<double, 2>(boost::extents[num_sinks][num_subsegments]));
    for (index iter = 0; iter < num_sinks; iter++)
    {
        for (index iter2 = 0; iter2 < num_subsegments; iter2++)
        {
            double distance = norm_2(mSubSegmentCoordinates[iter2].rGetLocation() - mSinkCoordinates[iter].rGetLocation());
            double term;
            if (distance <= equivalent_tissue_point_radius)
            {
                term = coefficient * (1.5 - 0.5 * (pow(distance / equivalent_tissue_point_radius, 2))) / equivalent_tissue_point_radius;
            }
            else
            {
                term = coefficient / distance;
            }
            (*p_interaction_matrix)[iter][iter2] = term;
        }
    }
    return p_interaction_matrix;
}

template<unsigned DIM>
boost::shared_ptr<boost::multi_array<double, 2> > GreensFunctionSolver<DIM>::GetVesselTissueInteractionMatrix()
{
    typedef boost::multi_array<double, 2>::index index;
    unsigned num_subsegments = mSubSegmentCoordinates.size();
    unsigned num_sinks = mSinkCoordinates.size();
    double coefficient = 1.0 / (4.0 * M_PI);

    double tissue_point_volume = pow(this->mpRegularGrid->GetSpacing(), 3);
    double equivalent_tissue_point_radius = pow(tissue_point_volume * 0.75 / M_PI, 0.333333);

    boost::shared_ptr<boost::multi_array<double, 2> > p_interaction_matrix(new boost::multi_array<double, 2>(boost::extents[num_subsegments][num_sinks]));
    for (index iter = 0; iter < num_subsegments; iter++)
    {
        for (index iter2 = 0; iter2 < num_sinks; iter2++)
        {
            double distance = norm_2(mSubSegmentCoordinates[iter].rGetLocation() - mSinkCoordinates[iter2].rGetLocation());
            double term;
            if (distance <= equivalent_tissue_point_radius)
            {
                term = coefficient * (1.5 - 0.5 * (pow(distance / equivalent_tissue_point_radius, 2)))/ equivalent_tissue_point_radius;
            }
            else
            {
                term = coefficient / distance;
            }
            (*p_interaction_matrix)[iter][iter2] = term;
        }
    }
    return p_interaction_matrix;
}

template<unsigned DIM>
void GreensFunctionSolver<DIM>::WriteSolution(std::map<std::string, std::vector<double> >& segmentPointData)
{
    // Write the tissue point data
    vtkSmartPointer<vtkXMLImageDataWriter> pImageDataWriter = vtkSmartPointer<vtkXMLImageDataWriter>::New();
    pImageDataWriter->SetFileName((this->mpOutputFileHandler->GetOutputDirectoryFullPath() + "/pde_solution.vti").c_str());
    pImageDataWriter->SetInput(this->mpVtkSolution);
    pImageDataWriter->Update();
    pImageDataWriter->Write();

    // Add the segment points
    vtkSmartPointer<vtkPolyData> pPolyData = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPoints> pPoints = vtkSmartPointer<vtkPoints>::New();
    for (unsigned i = 0; i < mSubSegmentCoordinates.size(); i++)
    {
        ChastePoint<DIM> location = mSubSegmentCoordinates[i];
        pPoints->InsertNextPoint(location[0], location[1], location[2]);
    }
    pPolyData->SetPoints(pPoints);

    // Add the segment point data
    std::map<std::string, std::vector<double> >::iterator segment_iter;
    for (segment_iter = segmentPointData.begin(); segment_iter != segmentPointData.end(); ++segment_iter)
    {
        vtkSmartPointer<vtkDoubleArray> pInfo = vtkSmartPointer<vtkDoubleArray>::New();
        pInfo->SetNumberOfComponents(1);
        pInfo->SetNumberOfTuples(mSubSegmentCoordinates.size());
        pInfo->SetName(segment_iter->first.c_str());

        for (unsigned i = 0; i < mSubSegmentCoordinates.size(); i++)
        {
            pInfo->SetValue(i, segmentPointData[segment_iter->first][i]);
        }
        pPolyData->GetPointData()->AddArray(pInfo);
    }

    vtkSmartPointer<vtkXMLPolyDataWriter> p_poldata_writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    p_poldata_writer->SetFileName((this->mpOutputFileHandler->GetOutputDirectoryFullPath() + "/segments.vtp").c_str());
    p_poldata_writer->SetInput(pPolyData);
    p_poldata_writer->Write();
}

// Explicit instantiation
template class GreensFunctionSolver<2>;
template class GreensFunctionSolver<3>;
