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

#include "VesselNetworkToMesh.hpp"
#include "NetworkToImage.hpp"
#include <vtkMarchingSquares.h>
#include <vtkCleanPolyData.h>
#include <vtkStripper.h>
#include <vtkSplineFilter.h>
#include <vtkTriangleFilter.h>
#include <vtkTransform.h>
#include <vtkClipPolyData.h>
#include <vtkPolyData.h>
#include <vtkIdList.h>
#include <vtkLine.h>
#include <vtkCellArray.h>
#include <vtkBox.h>
#include "ImageIO.hpp"
#include "VtkSurfaceWriter.hpp"
#include "UblasCustomFunctions.hpp"
#include "UblasIncludes.hpp"
#include "Debug.hpp"

template<unsigned DIM>
VesselNetworkToMesh<DIM>::VesselNetworkToMesh() :
    mpNetwork(),
    mMeshSize(10.0),
    mMeshDimension(3),
    mMesh2d(),
    mMesh3d()
{

}

template<unsigned DIM>
boost::shared_ptr<VesselNetworkToMesh<DIM> > VesselNetworkToMesh<DIM>::Create()
{
    MAKE_PTR(VesselNetworkToMesh<DIM>, pSelf);
    return pSelf;
}

template<unsigned DIM>
VesselNetworkToMesh<DIM>::~VesselNetworkToMesh()
{

}

template<unsigned DIM>
void VesselNetworkToMesh<DIM>::SetVesselNetwork(boost::shared_ptr<VesselNetwork<DIM> > pNetwork)
{
    mpNetwork = pNetwork;
}

template<unsigned DIM>
void VesselNetworkToMesh<DIM>::SetTargetMeshSize(double meshSize)
{
    mMeshSize = meshSize;
}

template<unsigned DIM>
void VesselNetworkToMesh<DIM>::SetMeshDimension(unsigned dimension)
{
    mMeshDimension = dimension;
}

template<unsigned DIM>
boost::shared_ptr<HybridMesh<DIM, DIM> > VesselNetworkToMesh<DIM>::Get2dMesh()
{
    return mMesh2d;
}

template<unsigned DIM>
void VesselNetworkToMesh<DIM>::DoMeshing2d()
{
    // Interpolate the network onto a regular grid
    boost::shared_ptr<NetworkToImage<DIM> > p_net_to_image = NetworkToImage<DIM>::Create();
    p_net_to_image->SetNetwork(mpNetwork);
    p_net_to_image->SetGridSpacing(2.0);
    p_net_to_image->SetPaddingFactors(0.1, 0.1, 0.0);
    p_net_to_image->SetImageDimension(2);
    p_net_to_image->Update();

    // Get the outer boundaries of the network
    vtkSmartPointer<vtkMarchingSquares> p_squares = vtkSmartPointer<vtkMarchingSquares>::New();
    p_squares->SetInput(p_net_to_image->GetOutput());
    p_squares->SetValue(0, 1);

    // Remove duplicate points and join up lines to form polylines
    vtkSmartPointer<vtkCleanPolyData> p_cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
    p_cleaner->SetInputConnection(p_squares->GetOutputPort());

    vtkSmartPointer<vtkStripper> p_stripper = vtkSmartPointer<vtkStripper>::New();
    p_stripper->SetInputConnection(p_cleaner->GetOutputPort());

    // Downsample and smooth the polyline
    vtkSmartPointer<vtkSplineFilter> p_spline= vtkSmartPointer<vtkSplineFilter>::New();
    p_spline->SetLength(10.0);
    p_spline->SetSubdivideToLength();
    p_spline->SetInputConnection(p_stripper->GetOutputPort());

    vtkSmartPointer<vtkTriangleFilter> p_triangle = vtkSmartPointer<vtkTriangleFilter>::New();
    p_triangle->SetInputConnection(p_spline->GetOutputPort());
    p_triangle->Update();

    vtkSmartPointer<vtkPolyData> p_cleaned = p_triangle->GetOutput();

    // Want flat ends on input and output network nodes. It is important to do this after smoothing
    std::vector<boost::shared_ptr<VesselNode<DIM> > > nodes = mpNetwork->GetNodes();
    c_vector<double, 3> box_axis = unit_vector<double>(3,0);

    for(unsigned idx=0; idx< nodes.size(); idx++)
    {
        if(nodes[idx]->GetFlowProperties()->IsInputNode() or nodes[idx]->GetFlowProperties()->IsOutputNode())
        {
            MARK;
            double radius = nodes[idx]->GetRadius();
            vtkSmartPointer<vtkBox> p_box = vtkSmartPointer<vtkBox>::New();
            p_box->SetBounds(-1.1*radius, 0.0, -1.1*radius, 1.1*radius, - 1.1*radius, 1.1*radius);

            c_vector<double, 3> loc;
            loc[0]= nodes[idx]->rGetLocation()[0];
            loc[1]= nodes[idx]->rGetLocation()[1];
            loc[2]=0.0;

            c_vector<double, 3> tangent;
            tangent[0]= nodes[idx]->GetSegment(0)->GetOppositeNode(nodes[idx])->rGetLocation()[0] - loc[0];
            tangent[1]= nodes[idx]->GetSegment(0)->GetOppositeNode(nodes[idx])->rGetLocation()[1] - loc[1];
            tangent[2] = 0.0;
            tangent /=norm_2(tangent);

            MARK;
            double rotation_angle = std::acos(inner_prod(box_axis, tangent))*(180.0/M_PI);
            c_vector<double, 3> rotation_axis = VectorProduct(box_axis, tangent);

            vtkSmartPointer<vtkTransform> p_tranform = vtkSmartPointer<vtkTransform>::New();
            if (std::abs(inner_prod(box_axis, tangent)) < 1.0 - 1.e-6)
            {
                p_tranform->RotateWXYZ(-rotation_angle, rotation_axis[0], rotation_axis[1], rotation_axis[2]);
            }
            else
            {
                p_tranform->RotateWXYZ(-rotation_angle, 0.0, 0.0, 1.0);
            }
            p_tranform->Translate(-loc[0],-loc[1], -loc[2]);
            p_box->SetTransform(p_tranform);

            MARK;
            vtkSmartPointer<vtkClipPolyData> p_clipper = vtkSmartPointer<vtkClipPolyData>::New();
            p_clipper->SetInput(p_cleaned);
            p_clipper->SetClipFunction(p_box);
            p_clipper->Update();

            VtkSurfaceWriter writer;
            writer.SetInput(p_clipper->GetOutput());
            writer.SetFileName("/home/grogan/test.vtp");
            writer.Write();

            MARK;
//            vtkSmartPointer<vtkPolyData> p_clipped;
//            p_clipped->DeepCopy(p_clipper->GetOutput());

            MARK;
            // Assuming we have two points with connectivity 1, join them
            std::vector<unsigned> edge_ids;
            vtkSmartPointer<vtkIdList> p_cell_list = vtkSmartPointer<vtkIdList>::New();
            for (unsigned idx =0; idx< p_clipper->GetOutput()->GetNumberOfPoints(); idx++)
            {
                p_clipper->GetOutput()->GetPointCells(idx, p_cell_list);
                if(p_cell_list->GetNumberOfIds()==1)
                {
                    edge_ids.push_back(idx);
                }
            }

            MARK;
            vtkSmartPointer<vtkLine> p_line = vtkSmartPointer<vtkLine>::New();
            p_line->GetPointIds()->SetId(0, edge_ids[0]);
            p_line->GetPointIds()->SetId(1, edge_ids[1]);
            p_clipper->GetOutput()->GetLines()->InsertNextCell(p_line);
            p_clipper->GetOutput()->Update();

            p_cleaned->DeepCopy(p_clipper->GetOutput());
            MARK;
        }
    }

    boost::shared_ptr<HybridMesh<DIM, DIM> > p_temp_mesh = HybridMesh<DIM, DIM>::Create();
    p_temp_mesh->GenerateTriMeshFromPolyData(p_cleaned);

    mMesh2d = p_temp_mesh;

//    tri_input_converter = chaste.mesh.converters.VtkToTriMesh()
//    tri_input_converter.input = cleaned
//    points, edges = tri_input_converter.update()
//
//    # First generate a coarse mesh and use it to probe for holes
//    mesh_info = MeshInfo()
//    mesh_info.set_points(points)
//    mesh_info.set_facets(edges)
//
//    # Allow for two different regions and holes
//    data = build(mesh_info)
//    points= data.points
//    edges = data.elements
//
//    vtk_points = vtk.vtkPoints()
//    for eachEdge in edges:
//        centx = (points[eachEdge[0]][0] + points[eachEdge[1]][0] + points[eachEdge[2]][0])/3.0
//        centy = (points[eachEdge[0]][1] + points[eachEdge[1]][1] + points[eachEdge[2]][1])/3.0
//        vtk_points.InsertNextPoint(centx, centy, 0.0)
//
//    # Get the values of the image data at the points
//    temp_polydata = vtk.vtkPolyData()
//    temp_polydata.SetPoints(vtk_points)
//    probe = vtk.vtkProbeFilter()
//    probe.SetInput(temp_polydata)
//    probe.SetSource(image)
//    probe.Update()
//    results = probe.GetOutput().GetPointData().GetScalars()
//
//    hole_locations = []
//    for idx in range(results.GetNumberOfTuples()):
//        if(results.GetTuple1(idx)==0):
//            loc = vtk_points.GetPoint(idx)
//            hole_locations.append([loc[0], loc[1]])
//
//    if(len(hole_locations)>0):
//        mesh_info.holes.resize(len(hole_locations))
//        for idx, eachHole in enumerate(hole_locations):
//            mesh_info.holes[idx] = [eachHole[0], eachHole[1]]
//
//    # Remesh with holes
//    data = build(mesh_info, max_volume = self.mesh_size)
//    tri_to_vtk = chaste.mesh.converters.TriMeshToVtkUnstructured()
//    tri_to_vtk.input = [data.points, data.elements]
//    mesh = tri_to_vtk.update()
//    return mesh
}

template<unsigned DIM>
void VesselNetworkToMesh<DIM>::DoMeshing3d()
{

}



template<unsigned DIM>
void VesselNetworkToMesh<DIM>::Update()
{
    if(!mpNetwork)
    {
        EXCEPTION("Vessel network required for meshing.");
    }

    if(mMeshDimension == 2)
    {
        DoMeshing2d();
    }
    else
    {
        DoMeshing3d();
    }
}

// Explicit instantiation
template class VesselNetworkToMesh<2> ;
template class VesselNetworkToMesh<3> ;
