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
//GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS unsignedERRUPTION)
//HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
//LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
//OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#include <stdio.h>
#include "Alarcon03HaematocritSolver.hpp"
#include "LinearSystem.hpp"
#include "VascularNode.hpp"
#include "Vessel.hpp"
#include "Exception.hpp"
#include "ReplicatableVector.hpp"

template<unsigned DIM>
Alarcon03HaematocritSolver<DIM>::Alarcon03HaematocritSolver(double haematocrit) : AbstractHaematocritSolver<DIM>(),
    mTHR(2.5),
    mAlpha(0.5),
    mHaematocrit(haematocrit)
{

}

template<unsigned DIM>
Alarcon03HaematocritSolver<DIM>::~Alarcon03HaematocritSolver()
{

}

template<unsigned DIM>
void Alarcon03HaematocritSolver<DIM>::SetTHR(double THR)
{
    mTHR = THR;
    assert(mTHR > 1);
}

template<unsigned DIM>
void Alarcon03HaematocritSolver<DIM>::SetAlpha(double Alpha)
{
    mAlpha = Alpha;
    assert(mAlpha < 1);
    assert(mAlpha > 0);
}

template<unsigned DIM>
void Alarcon03HaematocritSolver<DIM>::SetHaematocrit(double haematocrit)
{
    mHaematocrit = haematocrit;
}

template<unsigned DIM>
void Alarcon03HaematocritSolver<DIM>::Calculate(boost::shared_ptr<VascularNetwork<DIM> > vascularNetwork)
{

    // create extra data tables to aid with formation of coefficient matrix for haematocrit calculation
    std::vector<boost::shared_ptr<VascularNode<DIM> > > nodes = vascularNetwork->GetVesselEndNodes();
    std::vector< std::vector<unsigned> > VesselsFlowingOutOfNode(nodes.size());
    std::vector< std::vector<unsigned> > VesselsFlowingIntoNode(nodes.size());
    std::vector< std::vector<unsigned> > VesselsAttachedToNodeWithZeroFlow(nodes.size());

    for (unsigned i = 0; i < nodes.size(); i++)
    {
        for (unsigned j = 0; j < nodes[i]->GetNumberOfSegments(); j++)
        {
            boost::shared_ptr<Vessel<DIM> > p_vessel = nodes[i]->GetVesselSegment(j)->GetVessel();
            unsigned vessel_index = vascularNetwork->GetVesselIndex(p_vessel);
            double flow_rate = p_vessel->GetFlowRate();
            if (nodes[i] == p_vessel->GetStartNode())
            {
                if (flow_rate < 0)
                {
                    VesselsFlowingIntoNode[i].push_back(vessel_index);
                }
                else if (flow_rate > 0)
                {
                    VesselsFlowingOutOfNode[i].push_back(vessel_index);
                }
                else
                {
                    VesselsAttachedToNodeWithZeroFlow[i].push_back(vessel_index);
                }
            }
            else
            {
                if (flow_rate > 0)
                {
                    VesselsFlowingIntoNode[i].push_back(vessel_index);
                }
                else if (flow_rate < 0)
                {
                    VesselsFlowingOutOfNode[i].push_back(vessel_index);
                }
                else
                {
                    VesselsAttachedToNodeWithZeroFlow[i].push_back(vessel_index);
                }
            }
        }
    }

    // Set up the system
    unsigned number_of_vessels = vascularNetwork->GetNumberOfVessels();
    PetscInt lhsVectorSize = vascularNetwork->GetNumberOfVessels();
    unsigned pre_allocation_value;
    if(number_of_vessels < 3)
    {
        pre_allocation_value  = number_of_vessels;
    }
    else
    {
        pre_allocation_value = 3;
    }

    LinearSystem linearSystem(lhsVectorSize, pre_allocation_value);
    if(lhsVectorSize > 6)
    {
//        PetscOptionsSetValue("-pc_factor_mat_solver_package", "umfpack");
        linearSystem.SetPcType("hypre");
        linearSystem.SetKspType("preonly");
    }

    // Set the haematocrit of input vessels to the arterial level
    unsigned number_of_vessel_nodes = nodes.size();
    unsigned equation_number = 0;
    for (unsigned idx = 0; idx < number_of_vessel_nodes; idx++)
    {
        if (nodes[idx]->GetFlowProperties()->IsInputNode())
        {
            for (unsigned jdx = 0; jdx < nodes[idx]->GetNumberOfSegments(); jdx++)
            {
                linearSystem.AddToMatrixElement(equation_number, vascularNetwork->GetVesselIndex(nodes[idx]->GetVesselSegment(jdx)->GetVessel()), 1);
                linearSystem.SetRhsVectorElement(equation_number, mHaematocrit);
                equation_number++;
            }
        }
    }

    // Haematocrit conservation equations
    for (unsigned i = 0; i < number_of_vessel_nodes; i++)
    {
        unsigned number_of_inflow_vessels = VesselsFlowingIntoNode[i].size();
        unsigned number_of_outflow_vessels = VesselsFlowingOutOfNode[i].size();
        unsigned number_of_no_flow_vessels = VesselsAttachedToNodeWithZeroFlow[i].size();

        if (number_of_inflow_vessels + number_of_outflow_vessels + number_of_no_flow_vessels > 3)
        {
            EXCEPTION("The maximum number of coincident vessels at a node is 3.");
        }
        else if (number_of_inflow_vessels + number_of_outflow_vessels == 2)
        {
            if(!(number_of_inflow_vessels == 1 && number_of_outflow_vessels == 1))
            {
                EXCEPTION("Nodes with two conincident vessels must have one inflow and one outflow vessel.");
            }
            // output haematocrit equals input haematocrit
            linearSystem.AddToMatrixElement(equation_number, VesselsFlowingIntoNode[i][0], 1);
            linearSystem.AddToMatrixElement(equation_number, VesselsFlowingOutOfNode[i][0], -1);
            equation_number++;
        }
        else if (VesselsFlowingIntoNode[i].size() + VesselsFlowingOutOfNode[i].size() == 3)
        {
            if (VesselsFlowingIntoNode[i].size() == 2)
            {
                // output haematocrit equals sum of input haematocrits
                linearSystem.AddToMatrixElement(equation_number, VesselsFlowingIntoNode[i][0], 1);
                linearSystem.AddToMatrixElement(equation_number, VesselsFlowingIntoNode[i][1], 1);
                linearSystem.AddToMatrixElement(equation_number, VesselsFlowingOutOfNode[i][0], -1);
                equation_number++;
            }
            else if (VesselsFlowingIntoNode[i].size() == 1)
            {
                // bifurcation
                linearSystem.AddToMatrixElement(equation_number, VesselsFlowingIntoNode[i][0], 1);
                linearSystem.AddToMatrixElement(equation_number, VesselsFlowingOutOfNode[i][0], -1);
                linearSystem.AddToMatrixElement(equation_number, VesselsFlowingOutOfNode[i][1], -1);
                equation_number++;

                double radius0 = vascularNetwork->GetVessel(VesselsFlowingOutOfNode[i][0])->GetRadius();
                double radius1 = vascularNetwork->GetVessel(VesselsFlowingOutOfNode[i][1])->GetRadius();
                double out_flow_velocity0 = fabs(vascularNetwork->GetVessel(VesselsFlowingOutOfNode[i][0])->GetFlowRate())/(M_PI * radius0 * radius0);
                double out_flow_velocity1 = fabs(vascularNetwork->GetVessel(VesselsFlowingOutOfNode[i][1])->GetFlowRate())/(M_PI * radius1 * radius1);

                if (out_flow_velocity0 >= out_flow_velocity1)
                {
                    if (out_flow_velocity0 < mTHR*out_flow_velocity1)
                    {
                        linearSystem.AddToMatrixElement(equation_number, VesselsFlowingOutOfNode[i][0], 1);
                        linearSystem.AddToMatrixElement(equation_number, VesselsFlowingOutOfNode[i][1], -mAlpha*(out_flow_velocity0/out_flow_velocity1));
                    }
                    else
                    {
                        linearSystem.AddToMatrixElement(equation_number, VesselsFlowingOutOfNode[i][0], 1);
                        linearSystem.AddToMatrixElement(equation_number, VesselsFlowingIntoNode[i][0], -1);
                    }

                }
                else
                {
                    if (out_flow_velocity1 < mTHR*out_flow_velocity0)
                    {
                        linearSystem.AddToMatrixElement(equation_number, VesselsFlowingOutOfNode[i][1], 1);
                        linearSystem.AddToMatrixElement(equation_number, VesselsFlowingOutOfNode[i][0], -mAlpha*(out_flow_velocity1/out_flow_velocity0));
                    }
                    else
                    {
                        linearSystem.AddToMatrixElement(equation_number, VesselsFlowingOutOfNode[i][1], 1);
                        linearSystem.AddToMatrixElement(equation_number, VesselsFlowingIntoNode[i][0], -1);
                    }
                }
                equation_number++;
            }
        }
    }

    // zero flow vessels have zero haematocrit
    linearSystem.AssembleIntermediateLinearSystem();
    for (unsigned idx = 0; idx < number_of_vessels; idx++)
    {
        if (vascularNetwork->GetVessel(idx)->GetFlowRate() == 0)
        {
            linearSystem.AddToMatrixElement(equation_number, idx, 1.0);
            equation_number++;
        }
    }

    Vec solution = PetscTools::CreateVec(number_of_vessels);
     // Does an initial guess do anything with a direct solver?
    linearSystem.AssembleFinalLinearSystem();
    //linearSystem.DisplayMatrix();
    //linearSystem.DisplayRhs();
    solution = linearSystem.Solve();

    // deal with minor rounding errors in calculation
    ReplicatableVector a(solution);
    for (unsigned i = 0; i < number_of_vessels; i++)
    {
        if (a[i] < pow(10.0,-15))
        {
            a[i] = 0;
        }
    }

    // assign haematocrit levels to vessels
    for (unsigned i = 0; i < number_of_vessels; i++)
    {
        for (unsigned jdx = 0; jdx < vascularNetwork->GetVessel(i)->GetNumberOfSegments(); jdx++)
        {
            vascularNetwork->GetVessel(i)->GetSegment(jdx)->GetFlowProperties()->SetHaematocrit(double(a[i]));
        }
    }

    PetscTools::Destroy(solution);
}
// Explicit instantiation
template class Alarcon03HaematocritSolver<2>;
template class Alarcon03HaematocritSolver<3>;
