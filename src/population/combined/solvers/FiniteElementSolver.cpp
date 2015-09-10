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

#include <math.h>
#include "SimpleLinearEllipticSolver.hpp"
#include "VtkMeshWriter.hpp"
#include "ConstBoundaryCondition.hpp"
#include "FiniteElementSolver.hpp"
#include "VesselSurfaceGenerator.hpp"
#include "PlcMesh.hpp"
#include "CaVesselSegment.hpp"
#include "CaVascularNetwork.hpp"

template<unsigned DIM>
FiniteElementSolver<DIM>::FiniteElementSolver()
    : AbstractHybridSolver<DIM>(),
      mpDomain(),
      mGridSize(100.0),
      mpPde()
{

}

template<unsigned DIM>
FiniteElementSolver<DIM>::~FiniteElementSolver()
{

}

template<unsigned DIM>
void FiniteElementSolver<DIM>::SetDomain(boost::shared_ptr<Part<DIM> > pDomain)
{
    mpDomain = pDomain;
}

template<unsigned DIM>
void FiniteElementSolver<DIM>::SetMaxElementArea(double maxElementArea)
{
    mGridSize = maxElementArea;
}

template<unsigned DIM>
void FiniteElementSolver<DIM>::SetPde(boost::shared_ptr<HybridLinearEllipticPde<DIM, DIM> > pPde)
{
    mpPde = pPde;
}

template<unsigned DIM>
void FiniteElementSolver<DIM>::Solve(bool writeSolution)
{
    // If there is a vessel network add it to the domain
    if (this->mpNetwork)
    {
        mpDomain->AddVesselNetwork(this->mpNetwork, this->mBoundaryConditionType == BoundaryConditionType::SURFACE);
    }

    // Mesh the domain
    boost::shared_ptr<PlcMesh<DIM, DIM> > p_mesh = PlcMesh<DIM, DIM>::Create();
    p_mesh->GenerateFromPart(mpDomain, mGridSize);

    // Apply the boundary conditions
    double node_distance_tolerance = 1.e-3;
    BoundaryConditionsContainer<DIM, DIM, 1> bcc;

    if(!this->mBoundaryConditionType == BoundaryConditionType::SURFACE)
    {
        ConstBoundaryCondition<DIM>* p_fixed_boundary_condition = new ConstBoundaryCondition<DIM>(this->mBoundaryConditionValue);
        typename PlcMesh<DIM, DIM>::NodeIterator iter = p_mesh->GetNodeIteratorBegin();
        while (iter != p_mesh->GetNodeIteratorEnd())
        {
            c_vector<double, 3> location = zero_vector<double>(3);
            for(unsigned idx=0; idx<DIM; idx++)
            {
                location[idx] = (*iter).GetPoint()[idx];
            }
            std::pair<boost::shared_ptr<CaVesselSegment<DIM> >, double> seg_pair = this->mpNetwork->GetNearestSegment(location);
            double distance_to_segment = seg_pair.second;
            if (distance_to_segment < node_distance_tolerance)
            {
                bcc.AddDirichletBoundaryCondition(&(*iter), p_fixed_boundary_condition, 0, false);
            }
            ++iter;
        }
    }
    else
    {
        ConstBoundaryCondition<DIM>* p_fixed_boundary_condition = new ConstBoundaryCondition<DIM>(this->mBoundaryConditionValue);
        typename PlcMesh<DIM, DIM>::BoundaryNodeIterator iter = p_mesh->GetBoundaryNodeIteratorBegin();
        while (iter < p_mesh->GetBoundaryNodeIteratorEnd())
        {
            ChastePoint<DIM> location = (*iter)->GetPoint();
            std::pair<boost::shared_ptr<CaVesselSegment<DIM> >, double> seg_pair = this->mpNetwork->GetNearestSegment(location);
            if(seg_pair.second <= (seg_pair.first->GetRadius()+1.e-2))
            {
                bcc.AddDirichletBoundaryCondition(*iter, p_fixed_boundary_condition);
            }
            ++iter;
        }
    }

    // Do the solve
    SimpleLinearEllipticSolver<DIM, DIM> static_solver(p_mesh.get(), mpPde.get(), &bcc);
    ReplicatableVector static_solution_repl(static_solver.Solve());
    std::vector<double> output;
    for (unsigned idx = 0; idx < static_solution_repl.GetSize(); idx++)
    {
        output.push_back(static_solution_repl[idx]);
    }

    if(writeSolution)
    {
        this->Write(output, p_mesh);
    }

}

template<unsigned DIM>
void FiniteElementSolver<DIM>::Write(std::vector<double> output, boost::shared_ptr<PlcMesh<DIM, DIM> > p_mesh)
{
    // Write the output
    std::string fname;
    if(!this->mFilename.empty())
    {
        fname = this->mFilename;
    }
    else
    {
        fname = "solution";
    }

    VtkMeshWriter <DIM, DIM> mesh_writer(this->mWorkingDirectory, fname, false);
    if(output.size() > 0)
    {
        mesh_writer.AddPointData(mpPde->GetVariableName(), output);
    }
    mesh_writer.WriteFilesUsingMesh(*p_mesh);
}

// Explicit instantiation
template class FiniteElementSolver<2> ;
template class FiniteElementSolver<3> ;
