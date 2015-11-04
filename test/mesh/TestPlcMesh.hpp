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
 * Redistributions in binary form must reproduce the abovea copyright notice,
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

#ifndef TESTPLCMESH_HPP_
#define TESTPLCMESH_HPP_

#include <cxxtest/TestSuite.h>
#include <vector>
#include "SmartPointers.hpp"
#include "Polygon.hpp"
#include "Part.hpp"
#include "PlcMesh.hpp"
#include "VtkMeshWriter.hpp"
#include "OutputFileHandler.hpp"
#include "VasculatureGenerator.hpp"
#include "VascularNode.hpp"
#include "CaVesselSegment.hpp"
#include "CaVessel.hpp"
#include "CaVascularNetwork.hpp"

class TestPlcMesh : public CxxTest::TestSuite
{
private:

    boost::shared_ptr<CaVascularNetwork<3> > SetUpNetwork()
    {
        double vessel_length = 100;
        double radius = 10.0;
        double spacing = 3.0 * radius;
        unsigned num_vessels_per_row = 5;
        std::vector<boost::shared_ptr<VascularNode<3> > > start_nodes;
        std::vector<boost::shared_ptr<VascularNode<3> > > end_nodes;

        for(unsigned idx =0; idx<num_vessels_per_row; idx++)
        {
            for(unsigned jdx =0; jdx<num_vessels_per_row; jdx++)
            {
                double x_position = (spacing+2.0*radius) * double(idx) + spacing/2.0 + radius;
                double y_position = (spacing+2.0*radius) * double(jdx) + spacing/2.0 + radius;
                start_nodes.push_back(VascularNode<3>::Create(x_position, y_position, 0.0));
                end_nodes.push_back(VascularNode<3>::Create(x_position, y_position, vessel_length));
            }
        }

        std::vector<boost::shared_ptr<CaVessel<3> > > vessels;
        for(unsigned idx = 0; idx<start_nodes.size(); idx++)
        {
            start_nodes[idx]->SetRadius(radius);
            end_nodes[idx]->SetRadius(radius);
            vessels.push_back(CaVessel<3>::Create(CaVesselSegment<3>::Create(start_nodes[idx], end_nodes[idx])));
            vessels[idx]->GetSegment(0)->SetRadius(radius);
        }

        boost::shared_ptr<CaVascularNetwork<3> > p_network = CaVascularNetwork<3>::Create();
        p_network->AddVessels(vessels);
        return p_network;
    }

public:

    void dontTestMeshComposite2DPart() // Big memory problem, run at own risk!
    {
        boost::shared_ptr<Part<2> > p_part = Part<2>::Create();
        p_part->AddRectangle();
        p_part->AddCircle();

        boost::shared_ptr<PlcMesh<2> > p_mesh = PlcMesh<2>::Create();
        p_mesh->GenerateFromPart(p_part);

        VtkMeshWriter<2,2> mesh_writer("TestPlcMesh", "Composite2dMesh");
        mesh_writer.WriteFilesUsingMesh(*p_mesh);
    }

    void TestMeshCylinder()
    {
        boost::shared_ptr<Part<3> > p_part = Part<3>::Create();
        boost::shared_ptr<Polygon> p_circle = p_part->AddCircle();
        p_part->Extrude(p_circle);

        boost::shared_ptr<PlcMesh<3> > p_mesh = PlcMesh<3>::Create();
        p_mesh->GenerateFromPart(p_part, 20.0);

        VtkMeshWriter<3, 3> mesh_writer("TestPlcMesh", "Cylinder", false);
        mesh_writer.WriteFilesUsingMesh(*p_mesh);
    }

    void TestMeshCylinderWithVesselLine()
    {
        double vessel_length = 100.0;
        VasculatureGenerator<3> generator;
        boost::shared_ptr<CaVascularNetwork<3> > p_network = generator.GenerateSingleVessel(vessel_length);
        p_network->GetVessels()[0]->GetStartNode()->SetRadius(5.0);
        p_network->GetVessels()[0]->GetEndNode()->SetRadius(5.0);

        boost::shared_ptr<Part<3> > p_part = Part<3>::Create();
        boost::shared_ptr<Polygon> p_circle = p_part->AddCircle(100.0);
        p_part->Extrude(p_circle, 100.0);
        p_part->AddVesselNetwork(p_network);

        boost::shared_ptr<PlcMesh<3> > p_mesh = PlcMesh<3>::Create();
        p_mesh->GenerateFromPart(p_part, 100.0);

        VtkMeshWriter<3, 3> mesh_writer("TestPlcMesh", "CylinderWithVesselLine", false);
        mesh_writer.WriteFilesUsingMesh(*p_mesh);
    }

    void TestMeshCylinderWithVesselSurface()
    {
        double vessel_length = 100.0;
        VasculatureGenerator<3> generator;
        boost::shared_ptr<CaVascularNetwork<3> > p_network = generator.GenerateSingleVessel(vessel_length);
        p_network->GetVessels()[0]->GetStartNode()->SetRadius(5.0);
        p_network->GetVessels()[0]->GetEndNode()->SetRadius(5.0);

        boost::shared_ptr<Part<3> > p_part = Part<3>::Create();
        boost::shared_ptr<Polygon> p_circle = p_part->AddCircle(100.0);
        p_part->Extrude(p_circle, 100.0);
        p_part->AddVesselNetwork(p_network, true);

        boost::shared_ptr<PlcMesh<3> > p_mesh = PlcMesh<3>::Create();
        p_mesh->GenerateFromPart(p_part, 100.0);

        VtkMeshWriter < 3, 3 > mesh_writer("TestPlcMesh", "CylinderWithVesselSurface", false);
        mesh_writer.WriteFilesUsingMesh(*p_mesh);
    }

    void TestMeshCubeWithVesselLine()
    {
        double vessel_length = 100.0;
        VasculatureGenerator<3> generator;
        c_vector<double,3> centre = zero_vector<double>(3);
        centre[0] = vessel_length/2.0;
        centre[1] = vessel_length/2.0;
        boost::shared_ptr<CaVascularNetwork<3> > p_network = generator.GenerateSingleVessel(vessel_length, centre);
        p_network->GetVessels()[0]->GetStartNode()->SetRadius(5.0);
        p_network->GetVessels()[0]->GetEndNode()->SetRadius(5.0);

        boost::shared_ptr<Part<3> > p_part = Part<3>::Create();
        p_part->AddCuboid(vessel_length, vessel_length, vessel_length);
        p_part->AddVesselNetwork(p_network);

        boost::shared_ptr<PlcMesh<3> > p_mesh = PlcMesh<3>::Create();
        p_mesh->GenerateFromPart(p_part, 100.0);
        VtkMeshWriter <3, 3> mesh_writer("TestPlcMesh", "CubeWithVesselLine", false);
        mesh_writer.WriteFilesUsingMesh(*p_mesh);
    }

    void TestMeshCubeWithVesselSurface()
    {
        double vessel_length = 100.0;
        VasculatureGenerator<3> generator;
        c_vector<double,3> centre = zero_vector<double>(3);
        centre[0] = vessel_length/2.0;
        centre[1] = vessel_length/2.0;
        boost::shared_ptr<CaVascularNetwork<3> > p_network = generator.GenerateSingleVessel(vessel_length, centre);
        p_network->GetVessels()[0]->GetStartNode()->SetRadius(10.0);
        p_network->GetVessels()[0]->GetEndNode()->SetRadius(10.0);

        boost::shared_ptr<Part<3> > p_part = Part<3>::Create();
        p_part->AddCuboid(2.0 * vessel_length, 2.0 * vessel_length, vessel_length);
        p_part->AddVesselNetwork(p_network, true);

        boost::shared_ptr<PlcMesh<3> > p_mesh = PlcMesh<3>::Create();
        p_mesh->GenerateFromPart(p_part, 100.0);
        VtkMeshWriter<3, 3> mesh_writer("TestPlcMesh", "CubeWithVesselSurface", false);
        mesh_writer.WriteFilesUsingMesh(*p_mesh);
    }

    void TestMeshCubeWithVesselSurfaceInternal()
    {
        double vessel_length = 100.0;
        VasculatureGenerator<3> generator;
        c_vector<double,3> centre = zero_vector<double>(3);
        centre[0] = vessel_length/2.0;
        centre[1] = vessel_length/2.0;
        boost::shared_ptr<CaVascularNetwork<3> > p_network = generator.GenerateSingleVessel(vessel_length, centre);
        p_network->GetVessels()[0]->GetStartNode()->SetRadius(10.0);
        p_network->GetVessels()[0]->GetEndNode()->SetRadius(10.0);

        c_vector<double,3> translate = zero_vector<double>(3);
        translate[2] = -vessel_length/2.0;

        boost::shared_ptr<Part<3> > p_part = Part<3>::Create();
        p_part->AddCuboid(vessel_length, vessel_length, 2.0*vessel_length);
        p_part->Translate(translate);
        p_part->AddVesselNetwork(p_network, true);

        boost::shared_ptr<PlcMesh<3> > p_mesh = PlcMesh<3>::Create();
        p_mesh->GenerateFromPart(p_part, 100.0);
        VtkMeshWriter<3, 3> mesh_writer("TestPlcMesh", "CubeWithVesselSurface", false);
        mesh_writer.WriteFilesUsingMesh(*p_mesh);
    }

    void TestParrallelVesselSurfaceCube()
    {
        double vessel_length = 100;
        double radius = 10.0;
        double spacing = 3.0 * radius;
        unsigned num_vessels_per_row = 5;

        double domain_width = num_vessels_per_row * (spacing + 2.0* radius);
        double domain_height = num_vessels_per_row * (spacing + 2.0* radius);
        boost::shared_ptr<Part<3> > p_part = Part<3>::Create();
        p_part->AddCuboid(domain_width, domain_height, vessel_length);
        p_part->AddVesselNetwork(SetUpNetwork(), true);

        boost::shared_ptr<PlcMesh<3> > p_mesh = PlcMesh<3>::Create();
        p_mesh->GenerateFromPart(p_part, 1000.0);
        VtkMeshWriter < 3, 3 > mesh_writer("TestPlcMesh", "ParrallelVesselSurface", false);
        mesh_writer.WriteFilesUsingMesh(*p_mesh);
    }
};

#endif /*TESTPLCMESH_HPP_*/
