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

#include <sstream>
#include <boost/lexical_cast.hpp>
#ifdef CHASTE_VTK
#define _BACKWARD_BACKWARD_WARNING_H 1 //Cut out the strstream deprecated warning for now (gcc4.3)
#include <vtkDoubleArray.h>
#include <vtkCellData.h>
#include <vtkCellArray.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkLine.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkSmartPointer.h>
#endif // CHASTE_VTK
#include "SmartPointers.hpp"
#include "Exception.hpp"

#include "VasculatureGenerator.hpp"

template<unsigned DIM>
VasculatureGenerator<DIM>::VasculatureGenerator()
{
}

template<unsigned DIM>
VasculatureGenerator<DIM>::~VasculatureGenerator()
{
}

template<unsigned DIM>
void VasculatureGenerator<DIM>::PatternUnitByTranslation(boost::shared_ptr<CaVascularNetwork<DIM> > input_unit,
                                                         std::vector<unsigned> numberOfUnits)
{

    // Get unit dimensions
    std::vector<std::pair<double, double> > extents= input_unit->GetExtents();

    // For each number of units in each dimension copy the original vessels and move the copy to the desired location
    double num_x = 0;
    double num_y = 0;
    double num_z = 0;

    if(numberOfUnits.size() >= 1)
    {
        num_x = numberOfUnits[0];

        if(numberOfUnits.size() >= 2)
        {
            num_y = numberOfUnits[1];
            if(numberOfUnits.size() >= 3 && DIM ==3)
            {
                num_z = numberOfUnits[2];
            }
        }
    }

    // Keep track of the current vessels
    std::vector<boost::shared_ptr<CaVessel<DIM> > > original_vessels = input_unit->GetVessels();

    for(unsigned idx =0; idx < num_x; idx++)
    {
        c_vector<double, DIM> translation_vector;
        translation_vector[0] = double(idx+1) * (extents[0].second - extents[0].first);
        translation_vector[1] = 0.0;
        if(DIM==3)
        {
            translation_vector[2] = 0.0;
        }
        std::vector<boost::shared_ptr<CaVessel<DIM> > > copied_vessels = input_unit->CopyVessels(original_vessels);
        input_unit->Translate(translation_vector, copied_vessels);
    }

    input_unit->MergeCoincidentNodes();
    std::vector<boost::shared_ptr<CaVessel<DIM> > > x_transform_vessels = input_unit->GetVessels();

    for(unsigned idx =0; idx < num_y; idx++)
    {
        c_vector<double, DIM> translation_vector;
        translation_vector[0] = 0.0;
        translation_vector[1] = double(idx+1) * (extents[1].second - extents[1].first);
        if(DIM==3)
        {
            translation_vector[2] = 0.0;
        }
        std::vector<boost::shared_ptr<CaVessel<DIM> > > copied_vessels = input_unit->CopyVessels(x_transform_vessels);
        input_unit->Translate(translation_vector, copied_vessels);
    }
    input_unit->MergeCoincidentNodes();
    std::vector<boost::shared_ptr<CaVessel<DIM> > > y_transform_vessels = input_unit->GetVessels();

    for(unsigned idx =0; idx < num_z; idx++)
    {
        c_vector<double, DIM> translation_vector;
        translation_vector[0] = 0.0;
        translation_vector[1] = 0.0;
        if(DIM==3)
        {
            translation_vector[2] = double(idx+1) * (extents[2].second - extents[2].first);
        }
        std::vector<boost::shared_ptr<CaVessel<DIM> > > copied_vessels = input_unit->CopyVessels(y_transform_vessels);
        input_unit->Translate(translation_vector, copied_vessels);
    }
    input_unit->MergeCoincidentNodes();
}

template<unsigned DIM>
boost::shared_ptr<CaVascularNetwork<DIM> > VasculatureGenerator<DIM>::GenerateSimpleDivergeAndConvergeNetwork(c_vector<double, DIM> start_location,
                                                                                                              c_vector<double, DIM> end_location,
                                                                                                              double segmentLength,
                                                                                                              const std::string& rFileName)
{
    boost::shared_ptr<CaVascularNetwork<DIM> > p_network(new CaVascularNetwork<DIM>());

    // Create vessels at the start and end location
    std::vector<boost::shared_ptr<VascularNode<DIM> > > nodes;
    nodes.push_back(VascularNode<DIM>::Create(start_location));
    nodes[0]->SetLabel("Stalk");

    c_vector<double,DIM> network_direction = end_location - start_location;
    network_direction /= norm_2(network_direction);
    nodes.push_back(VascularNode<DIM>::Create(start_location + segmentLength*network_direction));
    nodes[1]->SetLabel("Tip");
    nodes.push_back(VascularNode<DIM>::Create(end_location));
    nodes[2]->SetLabel("Stalk");
    nodes.push_back(VascularNode<DIM>::Create(end_location - segmentLength*network_direction));
    nodes[3]->SetLabel("Tip");

    std::vector<boost::shared_ptr<CaVessel<DIM> > > start_vessels;
    start_vessels.push_back(CaVessel<DIM>::Create(CaVesselSegment<DIM>::Create(nodes[0], nodes[1])));
    start_vessels[0]->SetLabel("Start");

    std::vector<boost::shared_ptr<CaVessel<DIM> > > end_vessels;
    end_vessels.push_back(CaVessel<DIM>::Create(CaVesselSegment<DIM>::Create(nodes[2], nodes[3])));
    end_vessels[0]->SetLabel("End");

    p_network->AddVessels(start_vessels);
    p_network->AddVessels(end_vessels);

    // Start persistent random walks from both ends
    std::string tip_label = "Tip";
    for(unsigned outer_index=0; outer_index<10; outer_index++)
    {
        // Add segments to each vessel with a tip
        std::vector<boost::shared_ptr<CaVessel<DIM> > > vessels = p_network->GetVessels();
        for(unsigned idx=0; idx<vessels.size(); idx++)
        {
            // Get the start and end segments
            boost::shared_ptr<CaVesselSegment<DIM> > p_start_segment = vessels[idx]->GetSegment(0);
            boost::shared_ptr<CaVesselSegment<DIM> > p_end_segment = vessels[idx]->GetSegment(vessels[idx]->GetNumberOfSegments()-1);

            // Check if either have tips
            if(p_start_segment->GetNode(0)->rGetLabel().compare(tip_label) == 0)
            {
                // Get the segment vector
                c_vector<double, DIM> segment_vector = p_start_segment->GetNode(0)->GetLocationVector() - p_start_segment->GetNode(1)->GetLocationVector();
                segment_vector /= norm_2(segment_vector);

                // Add new segment
                boost::shared_ptr<VascularNode<DIM> > p_new_node = VascularNode<DIM>::Create(p_start_segment->GetNode(0)->GetLocationVector()
                                                                                             + segment_vector*segmentLength);
                p_new_node->SetLabel("Tip");
                p_start_segment->GetNode(0)->SetLabel("Stalk");
                vessels[idx]->AddSegment(CaVesselSegment<DIM>::Create(p_start_segment->GetNode(0),p_new_node));
            }

            else if(p_start_segment->GetNode(1)->rGetLabel().compare(tip_label) == 0)
            {
                // Get the segment vector
                c_vector<double, DIM> segment_vector = p_start_segment->GetNode(1)->GetLocationVector() - p_start_segment->GetNode(0)->GetLocationVector();
                segment_vector /= norm_2(segment_vector);

                // Add new segment
                boost::shared_ptr<VascularNode<DIM> > p_new_node = VascularNode<DIM>::Create(p_start_segment->GetNode(1)->GetLocationVector()
                                                                                             + segment_vector*segmentLength);
                p_new_node->SetLabel("Tip");
                p_start_segment->GetNode(1)->SetLabel("Stalk");
                vessels[idx]->AddSegment(CaVesselSegment<DIM>::Create(p_start_segment->GetNode(1),p_new_node));
            }

            if(p_end_segment->GetNode(0)->rGetLabel().compare(tip_label) == 0)
            {
                // Get the segment vector
                c_vector<double, DIM> segment_vector = p_end_segment->GetNode(0)->GetLocationVector() - p_end_segment->GetNode(1)->GetLocationVector();
                segment_vector /= norm_2(segment_vector);

                boost::shared_ptr<VascularNode<DIM> > p_new_node = VascularNode<DIM>::Create(p_end_segment->GetNode(0)->GetLocationVector()
                                                                                             + segment_vector*segmentLength);
                p_new_node->SetLabel("Tip");
                p_end_segment->GetNode(0)->SetLabel("Stalk");
                vessels[idx]->AddSegment(CaVesselSegment<DIM>::Create(p_end_segment->GetNode(0),p_new_node));
            }

            else if(p_end_segment->GetNode(1)->rGetLabel().compare(tip_label) == 0)
            {
                // Get the segment vector
                c_vector<double, DIM> segment_vector = p_end_segment->GetNode(1)->GetLocationVector() - p_end_segment->GetNode(0)->GetLocationVector();
                segment_vector /= norm_2(segment_vector);

                // Add new segment
                boost::shared_ptr<VascularNode<DIM> > p_new_node = VascularNode<DIM>::Create(p_end_segment->GetNode(1)->GetLocationVector()
                                                                                             + segment_vector*segmentLength);
                p_new_node->SetLabel("Tip");
                p_end_segment->GetNode(1)->SetLabel("Stalk");
                vessels[idx]->AddSegment(CaVesselSegment<DIM>::Create(p_end_segment->GetNode(1),p_new_node));
            }



        }
        unsigned cast_index = outer_index;
        std::string file_name = rFileName + "inc_" + boost::lexical_cast<std::string>(cast_index) + ".vtp";
        p_network->Write(file_name);
    }

    return p_network;
}

template<unsigned DIM>
boost::shared_ptr<CaVascularNetwork<DIM> > VasculatureGenerator<DIM>::GenerateHexagonalUnit(double vesselLength)
{
    // Generate the nodes
    std::vector<boost::shared_ptr<VascularNode<DIM> > > nodes;
    nodes.push_back(VascularNode<DIM>::Create(0.0, 0.0, 0.0)); //0
    nodes.push_back(VascularNode<DIM>::Create(vesselLength, vesselLength, 0.0)); //1
    nodes.push_back(VascularNode<DIM>::Create(0.0, 2.0 * vesselLength, 0.0)); //2
    nodes.push_back(VascularNode<DIM>::Create(2.0 * vesselLength, vesselLength, 0.0)); //3
    nodes.push_back(VascularNode<DIM>::Create(3.0 * vesselLength, 0.0, 0.0)); //4
    nodes.push_back(VascularNode<DIM>::Create(4.0 * vesselLength, 0.0, 0.0)); //5
    nodes.push_back(VascularNode<DIM>::Create(3.0 * vesselLength, 2.0 * vesselLength, 0.0)); //6
    if (DIM == 3)
    {
        nodes.push_back(VascularNode<DIM>::Create(0.0, 0.0, 1.5*vesselLength)); //7
        nodes.push_back(VascularNode<DIM>::Create(vesselLength, vesselLength, 1.5*vesselLength)); //8
        nodes.push_back(VascularNode<DIM>::Create(2.0 * vesselLength, vesselLength, 1.5*vesselLength)); //9
        nodes.push_back(VascularNode<DIM>::Create(3.0 * vesselLength, 0.0, 1.5*vesselLength)); //9
    }

    // Generate the segments and vessels
    std::vector<boost::shared_ptr<CaVessel<DIM> > > vessels;
    vessels.push_back(CaVessel<DIM>::Create(CaVesselSegment<DIM>::Create(nodes[0], nodes[1])));
    vessels.push_back(CaVessel<DIM>::Create(CaVesselSegment<DIM>::Create(nodes[1], nodes[2])));
    vessels.push_back(CaVessel<DIM>::Create(CaVesselSegment<DIM>::Create(nodes[1], nodes[3])));
    vessels.push_back(CaVessel<DIM>::Create(CaVesselSegment<DIM>::Create(nodes[3], nodes[4])));
    vessels.push_back(CaVessel<DIM>::Create(CaVesselSegment<DIM>::Create(nodes[4], nodes[5])));
    vessels.push_back(CaVessel<DIM>::Create(CaVesselSegment<DIM>::Create(nodes[3], nodes[6])));

    if (DIM == 3)
    {
        vessels.push_back(CaVessel<DIM>::Create(CaVesselSegment<DIM>::Create(nodes[0], nodes[7])));
        vessels.push_back(CaVessel<DIM>::Create(CaVesselSegment<DIM>::Create(nodes[1], nodes[8])));
        vessels.push_back(CaVessel<DIM>::Create(CaVesselSegment<DIM>::Create(nodes[3], nodes[9])));
        vessels.push_back(CaVessel<DIM>::Create(CaVesselSegment<DIM>::Create(nodes[4], nodes[10])));
    }

    // Generate the network
    boost::shared_ptr<CaVascularNetwork<DIM> > pNetwork(new CaVascularNetwork<DIM>());
    pNetwork->AddVessels(vessels);

    return pNetwork;
}

template<unsigned DIM>
boost::shared_ptr<CaVascularNetwork<DIM> > VasculatureGenerator<DIM>::GenerateBifurcationUnit(double vesselLength,
        c_vector<double, DIM> startPosition)
{
    // Generate the nodes
    std::vector<boost::shared_ptr<VascularNode<DIM> > > nodes;
    nodes.push_back(VascularNode<DIM>::Create(0.0, vesselLength, 0.0)); //0
    nodes.push_back(VascularNode<DIM>::Create(vesselLength, vesselLength, 0.0)); //1
    nodes.push_back(VascularNode<DIM>::Create(2.0 * vesselLength, 2.0 * vesselLength, 0.0)); //2
    nodes.push_back(VascularNode<DIM>::Create(2.0 * vesselLength, 0.0, 0.0)); //3
    nodes.push_back(VascularNode<DIM>::Create(3.0 * vesselLength, vesselLength, 0.0)); //4
    nodes.push_back(VascularNode<DIM>::Create(4.0 * vesselLength, vesselLength, 0.0)); //5

    // Generate the segments and vessels
    std::vector<boost::shared_ptr<CaVessel<DIM> > > vessels;
    vessels.push_back(CaVessel<DIM>::Create(CaVesselSegment<DIM>::Create(nodes[0], nodes[1])));
    vessels.push_back(CaVessel<DIM>::Create(CaVesselSegment<DIM>::Create(nodes[1], nodes[2])));
    vessels.push_back(CaVessel<DIM>::Create(CaVesselSegment<DIM>::Create(nodes[1], nodes[3])));
    vessels.push_back(CaVessel<DIM>::Create(CaVesselSegment<DIM>::Create(nodes[2], nodes[4])));
    vessels.push_back(CaVessel<DIM>::Create(CaVesselSegment<DIM>::Create(nodes[3], nodes[4])));
    vessels.push_back(CaVessel<DIM>::Create(CaVesselSegment<DIM>::Create(nodes[4], nodes[5])));

    // Generate the network
    boost::shared_ptr<CaVascularNetwork<DIM> > p_network(new CaVascularNetwork<DIM>());
    p_network->AddVessels(vessels);
    p_network->Translate(startPosition);
    return p_network;
}

template<unsigned DIM>
boost::shared_ptr<CaVascularNetwork<DIM> > VasculatureGenerator<DIM>::GenerateSingleVessel(double vesselLength,
        c_vector<double, DIM> startPosition)
{
    std::vector<boost::shared_ptr<VascularNode<DIM> > > nodes;
    nodes.push_back(VascularNode<DIM>::Create(0.0, 0.0, 0.0)); //0
    if(DIM==2)
    {
        nodes.push_back(VascularNode<DIM>::Create(0.0, vesselLength, 0.0)); //1
    }
    else
    {
        nodes.push_back(VascularNode<DIM>::Create(0.0, 0.0, vesselLength)); //1
    }

    std::vector<boost::shared_ptr<CaVessel<DIM> > > vessels;
    vessels.push_back(CaVessel<DIM>::Create(CaVesselSegment<DIM>::Create(nodes[0], nodes[1])));

    // Generate the network
    boost::shared_ptr<CaVascularNetwork<DIM> > p_network(new CaVascularNetwork<DIM>());
    p_network->AddVessels(vessels);
    p_network->Translate(startPosition);
    return p_network;
}

template<typename T>
std::string to_string(T const& value)
{
    std::stringstream sstr;
    sstr << value;
    return sstr.str();
}

template<unsigned DIM>
boost::shared_ptr<CaVascularNetwork<DIM> > VasculatureGenerator<DIM>::GenerateHexagonalNetwork(double width,
                                                                                               double height,
                                                                                               double vessel_length)
{

    // Determine the number of repeating hexagonal units in the x and y directions based on the input height, width and vessel length
    double horizontal_vessel_length = vessel_length;
    double diagonal_vessel_length = floor(vessel_length / (std::sqrt(2) * 2e-5) + 0.5) * 2e-5;
    unsigned units_in_x_direction = floor(
            ((width - horizontal_vessel_length - 2.0 * diagonal_vessel_length)
                    / (2.0 * (horizontal_vessel_length + diagonal_vessel_length))));
    unsigned units_in_y_direction = floor(((height) / (2.0 * diagonal_vessel_length)));

    // Ensure there are a minimal number of units required to generate a functional network
    if (units_in_x_direction <= 1 || units_in_y_direction <= 1)
    {
        std::string message =
                "Insufficient number of repeating units specified for the hexagonal network. Minimum length in x = ";
        message.append(
                to_string<double>(
                        2.0 * (horizontal_vessel_length + diagonal_vessel_length) + horizontal_vessel_length
                                + 2.0 * diagonal_vessel_length));
        message.append(". Minimum length in y = ");
        message.append(to_string<double>(2.0 * diagonal_vessel_length));
        message.append(".");
        EXCEPTION(message);
    }

    // Generate an array of vessels with no spatial info
    unsigned number_of_vessels = (units_in_x_direction * units_in_y_direction * 6) + (units_in_y_direction * 5)
            + units_in_x_direction;
    boost::shared_ptr<CaVascularNetwork<DIM> > pVesselNetwork(new CaVascularNetwork<DIM>());
    std::vector<boost::shared_ptr<CaVessel<DIM> > > vessel_array;
    std::vector<boost::shared_ptr<CaVesselSegment<DIM> > > vessel_segment_array;

    for (unsigned i = 0; i < number_of_vessels; i++)
    {
        // instantiate set of vessel segments with any nodes then later change location of that node
        // with GetNode(index)->SetLocation(ChastePoint ... )
        boost::shared_ptr<VascularNode<DIM> > node1(VascularNode<DIM>::Create(0, 0, 0));
        boost::shared_ptr<VascularNode<DIM> > node2(VascularNode<DIM>::Create(0, 0, 0));
        vessel_segment_array.push_back(CaVesselSegment<DIM>::Create(node1, node2));
    }

    // Set the left side coordinates of vessels in the network
    vessel_segment_array[0]->GetNode(0)->SetLocation(ChastePoint<DIM>(0.0, 0.0));
    vessel_segment_array[1]->GetNode(0)->SetLocation(
            ChastePoint<DIM>(diagonal_vessel_length + horizontal_vessel_length, diagonal_vessel_length));
    vessel_segment_array[2]->GetNode(0)->SetLocation(
            ChastePoint<DIM>(2.0 * diagonal_vessel_length + horizontal_vessel_length, 0.0));
    for (unsigned i = 3; i < (2 + 3 * units_in_x_direction); i++)
    {
        vessel_segment_array[i]->GetNode(0)->SetLocation(
                ChastePoint<DIM>(
                        vessel_segment_array[i - 3]->GetNode(0)->GetLocation()[0]
                                + 2.0 * (diagonal_vessel_length + horizontal_vessel_length),
                        vessel_segment_array[i - 3]->GetNode(0)->GetLocation()[1]));
    }
    vessel_segment_array[2 + 3 * units_in_x_direction]->GetNode(0)->SetLocation(
            ChastePoint<DIM>(0, 2.0 * diagonal_vessel_length));
    vessel_segment_array[2 + 3 * units_in_x_direction + 1]->GetNode(0)->SetLocation(
            ChastePoint<DIM>(diagonal_vessel_length, diagonal_vessel_length));
    vessel_segment_array[2 + 3 * units_in_x_direction + 2]->GetNode(0)->SetLocation(
            ChastePoint<DIM>(diagonal_vessel_length + horizontal_vessel_length, diagonal_vessel_length));

    for (unsigned i = (2 + 3 * units_in_x_direction + 3);
            i < (2 + 3 * units_in_x_direction + 3 * (units_in_x_direction + 1)); i++)
    {
        vessel_segment_array[i]->GetNode(0)->SetLocation(
                ChastePoint<DIM>(
                        vessel_segment_array[i - 3]->GetNode(0)->GetLocation()[0]
                                + 2.0 * (diagonal_vessel_length + horizontal_vessel_length),
                        vessel_segment_array[i - 3]->GetNode(0)->GetLocation()[1]));
    }

    for (unsigned i = (2 + 3 * units_in_x_direction + 3 * (units_in_x_direction + 1));
            i < number_of_vessels - units_in_x_direction; i++)
    {
        vessel_segment_array[i]->GetNode(0)->SetLocation(
                ChastePoint<DIM>(
                        vessel_segment_array[i - (2 + 3 * units_in_x_direction + 3 * (units_in_x_direction + 1))]->GetNode(
                                0)->GetLocation()[0],
                        vessel_segment_array[i - (2 + 3 * units_in_x_direction + 3 * (units_in_x_direction + 1))]->GetNode(
                                0)->GetLocation()[1] + 2.0 * diagonal_vessel_length));
    }
    vessel_segment_array[number_of_vessels - units_in_x_direction]->GetNode(0)->SetLocation(
            ChastePoint<DIM>(
                    2.0 * diagonal_vessel_length + horizontal_vessel_length,
                    vessel_segment_array[number_of_vessels - units_in_x_direction - 1]->GetNode(0)->GetLocation()[1]
                            + diagonal_vessel_length));

    for (unsigned i = 1; i < units_in_x_direction; i++)
    {
        vessel_segment_array[number_of_vessels - units_in_x_direction + i]->GetNode(0)->SetLocation(
                ChastePoint<DIM>(
                        vessel_segment_array[number_of_vessels - units_in_x_direction + i - 1]->GetNode(0)->GetLocation()[0]
                                + 2 * (diagonal_vessel_length + horizontal_vessel_length),
                        vessel_segment_array[number_of_vessels - units_in_x_direction + i - 1]->GetNode(0)->GetLocation()[1]));
    }

    // Set the right side coordinates of vessels in the network
    vessel_segment_array[0]->GetNode(1)->SetLocation(ChastePoint<DIM>(diagonal_vessel_length, diagonal_vessel_length));
    vessel_segment_array[1]->GetNode(1)->SetLocation(
            ChastePoint<DIM>(2.0 * diagonal_vessel_length + horizontal_vessel_length, 0.0));
    vessel_segment_array[2]->GetNode(1)->SetLocation(
            ChastePoint<DIM>(2.0 * (diagonal_vessel_length + horizontal_vessel_length), 0.0));

    for (unsigned i = 3; i < (2 + 3 * units_in_x_direction); i++)
    {
        vessel_segment_array[i]->GetNode(1)->SetLocation(
                ChastePoint<DIM>(
                        vessel_segment_array[i - 3]->GetNode(1)->GetLocation()[0]
                                + 2 * (diagonal_vessel_length + horizontal_vessel_length),
                        vessel_segment_array[i - 3]->GetNode(1)->GetLocation()[1]));
    }
    vessel_segment_array[2 + 3 * units_in_x_direction]->GetNode(1)->SetLocation(
            ChastePoint<DIM>(diagonal_vessel_length, diagonal_vessel_length));
    vessel_segment_array[2 + 3 * units_in_x_direction + 1]->GetNode(1)->SetLocation(
            ChastePoint<DIM>(diagonal_vessel_length + horizontal_vessel_length, diagonal_vessel_length));
    vessel_segment_array[2 + 3 * units_in_x_direction + 2]->GetNode(1)->SetLocation(
            ChastePoint<DIM>(2 * diagonal_vessel_length + horizontal_vessel_length, 2 * diagonal_vessel_length));

    for (unsigned i = (2 + 3 * units_in_x_direction + 3);
            i < (2 + 3 * units_in_x_direction + 3 * (units_in_x_direction + 1)); i++)
    {
        vessel_segment_array[i]->GetNode(1)->SetLocation(
                ChastePoint<DIM>(
                        vessel_segment_array[i - 3]->GetNode(1)->GetLocation()[0]
                                + 2 * (diagonal_vessel_length + horizontal_vessel_length),
                        vessel_segment_array[i - 3]->GetNode(1)->GetLocation()[1]));
    }

    for (unsigned i = (2 + 3 * units_in_x_direction + 3 * (units_in_x_direction + 1));
            i < number_of_vessels - units_in_x_direction; i++)
    {
        vessel_segment_array[i]->GetNode(1)->SetLocation(
                ChastePoint<DIM>(
                        vessel_segment_array[i - (2 + 3 * units_in_x_direction + 3 * (units_in_x_direction + 1))]->GetNode(
                                1)->GetLocation()[0],
                        vessel_segment_array[i - (2 + 3 * units_in_x_direction + 3 * (units_in_x_direction + 1))]->GetNode(
                                1)->GetLocation()[1] + 2 * diagonal_vessel_length));
    }

    vessel_segment_array[number_of_vessels - units_in_x_direction]->GetNode(1)->SetLocation(
            ChastePoint<DIM>(
                    2 * (diagonal_vessel_length + horizontal_vessel_length),
                    vessel_segment_array[number_of_vessels - units_in_x_direction - 1]->GetNode(0)->GetLocation()[1]
                            + diagonal_vessel_length));

    for (unsigned i = 1; i < units_in_x_direction; i++)
    {
        vessel_segment_array[number_of_vessels - units_in_x_direction + i]->GetNode(1)->SetLocation(
                ChastePoint<DIM>(
                        vessel_segment_array[number_of_vessels - units_in_x_direction + i - 1]->GetNode(1)->GetLocation()[0]
                                + 2 * (diagonal_vessel_length + horizontal_vessel_length),
                        vessel_segment_array[number_of_vessels - units_in_x_direction + i - 1]->GetNode(1)->GetLocation()[1]));
    }

    typename std::vector<boost::shared_ptr<CaVesselSegment<DIM> > >::iterator segment_iterator;

    for (unsigned i = 0; i < vessel_segment_array.size(); i++)
    {
        vessel_array.push_back(CaVessel<DIM>::Create(vessel_segment_array[i]));
    }

    pVesselNetwork->AddVessels(vessel_array);

    pVesselNetwork->MergeCoincidentNodes();

    return pVesselNetwork;
}

#ifdef CHASTE_VTK
template<unsigned DIM>
boost::shared_ptr<CaVascularNetwork<DIM> > VasculatureGenerator<DIM>::GenerateNetworkFromVtkFile(
        const std::string& filename)
{
    // Create an empty vessel network
    boost::shared_ptr<CaVascularNetwork<DIM> > pVesselNetwork(new CaVascularNetwork<DIM>());

    // Create a VTK PolyData object based on the contents of the input VTK file
    vtkSmartPointer<vtkXMLPolyDataReader> reader = vtkSmartPointer<vtkXMLPolyDataReader> ::New();
    vtkSmartPointer<vtkPolyData> pPolyData = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPointData> pPointData = vtkSmartPointer<vtkPointData>::New();

    reader->SetFileName(filename.c_str());
    reader->Update();
    pPolyData = reader->GetOutput();

    std::vector<double> radii;

    // Create the nodes
    std::vector<boost::shared_ptr<VascularNode<DIM> > > nodes;
    for (vtkIdType i = 0; i < pPolyData->GetNumberOfPoints(); i++)
    {
        double point_coords[3];
        pPolyData->GetPoint(i, point_coords);
        if (DIM < 3)
        {
            nodes.push_back(VascularNode<DIM>::Create(point_coords[0], point_coords[1]));
        }
        else
        {
            nodes.push_back(VascularNode<DIM>::Create(point_coords[0], point_coords[1], point_coords[2]));
        }
    }

    // Extract radii corresponding to each node from the VTK Polydata and store them in a list.
    pPointData = pPolyData->GetPointData();
    std::string radius_label1 = "Node Radius";
    std::string radius_label2 = "Radius";
    std::string radius_label3 = "radius";
    for (vtkIdType i = 0; i < pPointData->GetNumberOfArrays(); i++)
    {
        std::string array_name = pPointData->GetArrayName(i);
        if (array_name.compare(radius_label1) == 0 || array_name.compare(radius_label2) == 0||array_name.compare(radius_label3) == 0)
        {
            vtkDoubleArray* pScalars = vtkDoubleArray::SafeDownCast(pPointData->GetArray(i));
            for (vtkIdType j = 0; j < pScalars->GetNumberOfTuples(); j++)
            {
                radii.push_back(pScalars->GetValue(j));
            }
        }
    }

    // Extract vessels from the VTK Polydata. This is done by iterating over a VTK CellArray object which
    // returns a 'pointList' vtkIdList object. This object contains the point IDs of the nodes which make up
    // the vessel.
    vtkSmartPointer<vtkCellArray> pCellArray = vtkSmartPointer<vtkCellArray>::New();
    pCellArray = pPolyData->GetLines();

    ///\ todo Deleting this pointer at the end of the block caused problems, not sure why. Should check best way to delete
    // vtk pointers. Smart pointers don't seem to work for vtkIdType in v5.2.
    vtkIdType* pSegmentList;
    vtkIdType num_segments;
    for (int i = 0; i < pPolyData->GetNumberOfLines(); i++)
    {
        // Make a new vessel
        pCellArray->GetNextCell(num_segments, pSegmentList);
        std::vector<boost::shared_ptr<CaVesselSegment<DIM> > > segments;
        // Add segments to the vessels in order
        for (int j = 1; j < num_segments; j++)
        {
            boost::shared_ptr<CaVesselSegment<DIM> > p_segment = CaVesselSegment<DIM>::Create(nodes[pSegmentList[j - 1]],nodes[pSegmentList[j]]);
            if(unsigned(radii.size())>= pSegmentList[j])
            {
                p_segment->SetRadius(radii[pSegmentList[j]]);
            }
            segments.push_back(p_segment);
        }
        // Add the resulting vessel to the network
        pVesselNetwork->AddVessel(CaVessel<DIM>::Create(segments));
    }
    return pVesselNetwork;
}
#endif // CHASTE_VTK

//Explicit instantiation
template class VasculatureGenerator<2> ;
template class VasculatureGenerator<3> ;
