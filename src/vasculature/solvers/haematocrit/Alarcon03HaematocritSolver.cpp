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
// */
//
//#include "Alarcon03HaematocritCalculator.hpp"
//
//template<unsigned DIM>
//Alarcon03HaematocritCalculator<DIM>::Alarcon03HaematocritCalculator(double haematocrit)
//	: AbstractHaematocritCalculator<DIM>(),
//	  mTHR(2.5),
//	  mAlpha(0.5),
//	  mHaematocrit(haematocrit)
//{
//
//}
//
//template<unsigned DIM>
//Alarcon03HaematocritCalculator<DIM>::~Alarcon03HaematocritCalculator()
//{
//
//}
//
//template<unsigned DIM>
//void Alarcon03HaematocritCalculator<DIM>::SetTHR(double thr)
//{
//    mTHR = thr;
//    assert(mTHR > 1);
//}
//
//template<unsigned DIM>
//void Alarcon03HaematocritCalculator<DIM>::SetAlpha(double alpha)
//{
//    mAlpha = alpha;
//    assert(mAlpha < 1);
//    assert(mAlpha > 0);
//}
//
//void ConstantHaematocritSolver<DIM>::SetHaematocrit(double haematocrit)
//{
//	mHaematocrit = haematocrit;
//}
//
//template<unsigned DIM>
//void Alarcon03HaematocritCalculator<DIM>::Calculate(boost::shared_ptr<CaVascularNetwork<DIM> > vascularNetwork)
//{
//
//    // create extra data tables to aid with formation of coefficient matrix for haematocrit Calculator
//	std::vector<boost::shared_ptr<VascularNode<DIM> > > nodes = vascularNetwork->GetVectorOfNodes();
//	std::vector<boost::shared_ptr<CaVesselSegment<DIM> > > segments = vascularNetwork->GetVesselSegments();
//	unsigned num_nodes = nodes.size();
//
//    std::vector< std::vector<unsigned> > VesselsFlowingOutOfNode(num_nodes);
//    std::vector< std::vector<unsigned> > VesselsFlowingInToNode(num_nodes);
//    std::vector< std::vector<unsigned> > VesselsAttachedToNodeWithZeroFlow(num_nodes);
//
//    for (int node_index = 0; node_index < num_nodes; node_index++)
//    {
//    	boost::shared_ptr<VascularNode<DIM> > p_each_node = nodes[node_index];
//    	unsigned num_segments_on_node = p_each_node->GetNumberOfSegments();
//
//        for (int segment_index = 0; segment_index < num_segments_on_node; segment_index++)
//        {
//			boost::shared_ptr<CaVesselSegment<DIM> > p_each_segment = p_each_node->GetVesselSegment(segment_index);
//			typename std::vector<boost::shared_ptr<CaVesselSegment<DIM> > >::iterator segment_iter = std::find(segments.begin(), segments.end(), p_each_segment);
//			unsigned segment_global_index = std::distance(segments.begin(), segment_iter);
//
//			double flow_rate = p_each_segment->template GetData<double>("Flow Rate");
//
//			// Get the segment start and end nodes
//			boost::shared_ptr<VascularNode<DIM> > p_segment_start_node = p_each_segment->GetNode(0);
//			boost::shared_ptr<VascularNode<DIM> > p_segment_end_node = p_each_segment->GetNode(1);
//
//
//            if (flow_rate < 0.0)
//            {
//            	if (p_each_node == p_segment_start_node)
//            	{
//            		VesselsFlowingInToNode[node_index].push_back(segment_global_index);
//            	}
//            	else
//            	{
//            		VesselsFlowingOutOfNode[node_index].push_back(segment_global_index);
//            	}
//            }
//            else if (flow_rate > 0.0)
//            {
//            	if (p_each_node == p_segment_start_node)
//            	{
//            		VesselsFlowingOutOfNode[node_index].push_back(segment_global_index);
//            	}
//            	else
//            	{
//            		VesselsFlowingInToNode[node_index].push_back(segment_global_index);
//            	}
//            }
//            else
//            {
//            	VesselsAttachedToNodeWithZeroFlow[node_index].push_back(segment_global_index);
//            }
//        }
//    }
//
//    my::matrix CoeffMat(num_nodes, num_nodes);
//    int EquationNumber = 0;
//
//    // set up b vector for Calculator
//    vectordouble Haematocrit_bVector(num_nodes);
//
//    // equations which say that arterial input vessels have an arterial haematocrit level
//    for (unsigned node_index = 0; node_index < num_nodes; node_index++)
//    {
//        if (nodes[node_index]->template GetData<bool>("Is Input"))
//        {
//			boost::shared_ptr<CaVesselSegment<DIM> > p_each_segment = p_each_node->GetVesselSegment(0);
//			typename std::vector<boost::shared_ptr<CaVesselSegment<DIM> > >::iterator segment_iter = std::find(segments.begin(), segments.end(), p_each_segment);
//			unsigned segment_global_index = std::distance(segments.begin(), segment_iter);
//
//            CoeffMat.x[EquationNumber][segment_global_index] = 1;
//            Haematocrit_bVector.x[EquationNumber] = mHaematocrit;
//            EquationNumber++;
//        }
//    }
//
//    for (unsigned node_index = 0; node_index < num_nodes; node_index++)
//    {
//        unsigned num_inflow_vessels = VesselsFlowingInToNode[node_index].size();
//        unsigned num_outflow_vessels = VesselsFlowingOutOfNode[node_index].size();
//        unsigned num_no_flow_vessels = VesselsAttachedToNodeWithZeroFlow[node_index].size();
//
//        unsigned num_attached_vessels = num_inflow_vessels + num_outflow_vessels + num_no_flow_vessels;
//        if (num_attached_vessels > 3)
//        {
//        	EXCEPTION("No more than 3 vessels can intersect at a node.");
//        }
//        else if (num_inflow_vessels + num_outflow_vessels == 2)
//        {
//            // must be one flow going in to node and one flowing out for conservation
//
//        	if (num_inflow_vessels !=1 || num_outflow_vessels !=1)
//        	{
//        		EXCEPTION("Must have one in-flow and one out-flow vessel.");
//        	}
//
//            CoeffMat.x[EquationNumber][VesselsFlowingInToNode[node_index][0]] = 1;
//            CoeffMat.x[EquationNumber][VesselsFlowingOutOfNode[node_index][0]] = -1;
//            EquationNumber++;
//        }
//        else if (num_inflow_vessels + num_outflow_vessels == 3)
//        {
//
//            if (num_inflow_vessels == 2)
//            {
//                // conservation equation for node
//                CoeffMat.x[EquationNumber][VesselsFlowingInToNode[node_index][0]] = 1;
//                CoeffMat.x[EquationNumber][VesselsFlowingInToNode[node_index][1]] = 1;
//                CoeffMat.x[EquationNumber][VesselsFlowingOutOfNode[node_index][0]] = -1;
//                EquationNumber++;
//            }
//            if (num_inflow_vessels == 1)
//            {
//                // conservation equation for node
//                CoeffMat.x[EquationNumber][VesselsFlowingInToNode[node_index][0]] = 1;
//                CoeffMat.x[EquationNumber][VesselsFlowingOutOfNode[node_index][0]] = -1;
//                CoeffMat.x[EquationNumber][VesselsFlowingOutOfNode[node_index][1]] = -1;
//                EquationNumber++;
//
//                double absolute_flow_velocity_0 = segments[VesselsFlowingOutOfNode[node_index][0]]->template GetData<double>("Absolute Flow Velocity");
//                double absolute_flow_velocity_1 = segments[VesselsFlowingOutOfNode[node_index][1]]->template GetData<double>("Absolute Flow Velocity");
//
//                if (absolute_flow_velocity_0 >= absolute_flow_velocity_1)
//                {
//                    if (absolute_flow_velocity_0 < mTHR*absolute_flow_velocity_1)
//                    {
//                        CoeffMat.x[EquationNumber][VesselsFlowingOutOfNode[node_index][0]] = 1;
//                        CoeffMat.x[EquationNumber][VesselsFlowingOutOfNode[node_index][1]] = -mAlpha*(absolute_flow_velocity_0/absolute_flow_velocity_1);
//                    }
//                    else
//                    {
//                        CoeffMat.x[EquationNumber][VesselsFlowingOutOfNode[node_index][0]] = 1;
//                        CoeffMat.x[EquationNumber][VesselsFlowingInToNode[node_index][0]] = -1;
//                    }
//                }
//
//                else
//                {
//                    if (absolute_flow_velocity_1 < mTHR*absolute_flow_velocity_0)
//                    {
//                        CoeffMat.x[EquationNumber][VesselsFlowingOutOfNode[node_index][1]] = 1;
//                        CoeffMat.x[EquationNumber][VesselsFlowingOutOfNode[node_index][0]] = -mAlpha*(absolute_flow_velocity_1/absolute_flow_velocity_0);
//                    }
//                    else
//                    {
//                        CoeffMat.x[EquationNumber][VesselsFlowingOutOfNode[node_index][1]] = 1;
//                        CoeffMat.x[EquationNumber][VesselsFlowingInToNode[node_index][0]] = -1;
//                    }
//                }
//                EquationNumber++;
//            }
//        }
//    }
//
//	for (unsigned segment_index = 0; segment_index < segments.size(); segment_index++)
//	{
//		if (segments[segment_index]->template GetData<double>("Flow Velocity") == 0.0)
//		{
//            CoeffMat.x[EquationNumber][segment_index] = 1;
//            EquationNumber++;
//		}
//	}
//
//    // perform Calculator
//    vectordouble HaematocritLevels(segments.size());
//    HaematocritLevels = Haematocrit_bVector/CoeffMat;
//
//    // deal with minor rounding errors in Calculator
//    for (unsigned segment_index = 0; segment_index < segments.size(); segment_index++)
//    {
//        if (HaematocritLevels.x[segment_index] < pow(10.0,-15))
//        {
//            HaematocritLevels.x[segment_index] = 0.0;
//        }
//    }
//
//    // assign haematocrit levels to vessels
//    for (unsigned segment_index = 0; segment_index < segments.size(); segment_index++)
//    {
//    	segments[segment_index]->SetData["Haematocrit", HaematocritLevels.x[node_index]];
//    }
//}
//
//// Explicit instantiation
//
//template class Alarcon03HaematocritCalculator<2>;
//template class Alarcon03HaematocritCalculator<3>;
