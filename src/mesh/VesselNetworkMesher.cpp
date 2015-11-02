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

#define _BACKWARD_BACKWARD_WARNING_H 1 //Cut out the vtk deprecated warning for now (gcc4.3)
#include <vtkBox.h>
#include <vtkCylinder.h>
#include <vtkImplicitBoolean.h>
#include <vtkImplicitFunction.h>
#include <vtkSampleFunction.h>
#include <vtkContourFilter.h>
#include <vtkTransform.h>
#include <vtkBox.h>
#include <vtkSphere.h>
#include <stdlib.h>
#include "Exception.hpp"
#include "CaVessel.hpp"
#include "VascularNode.hpp"
#include "CaVesselSegment.hpp"
#include "UblasCustomFunctions.hpp"
#include "UblasIncludes.hpp"
#include "Vertex.hpp"
#include "Part.hpp"
#include "GeometryTools.hpp"

#include "VesselNetworkMesher.hpp"

template<unsigned DIM>
VesselNetworkMesher<DIM>::VesselNetworkMesher(boost::shared_ptr<CaVascularNetwork<DIM> > pVesselNetwork) :
        mpVesselNetwork(pVesselNetwork),
        mpSurface(vtkSmartPointer<vtkPolyData>::New())
{
}

template<unsigned DIM>
VesselNetworkMesher<DIM>::~VesselNetworkMesher()
{
}

template<unsigned DIM>
vtkSmartPointer<vtkPolyData> VesselNetworkMesher<DIM>::GetVtkSurface()
{
    c_vector<double, DIM> y_axis = unit_vector<double>(DIM,1);

    // Construct each segment as truncated cylinder or cone with spheres on the ends.
    std::vector<boost::shared_ptr<CaVesselSegment<DIM> > > segments = mpVesselNetwork->GetVesselSegments();

    std::vector<vtkSmartPointer<vtkImplicitFunction> > segment_functions;
    for(unsigned idx=0; idx<segments.size(); idx++)
    {
        // Check if the nodal radii are the same
        double start_radius = segments[idx]->GetNode(0)->GetRadius();
        double end_radius = segments[idx]->GetNode(1)->GetRadius();
        c_vector<double, DIM> start_location = segments[idx]->GetNode(0)->GetLocationVector();
        c_vector<double, DIM> end_location = segments[idx]->GetNode(1)->GetLocationVector();

        if(std::abs(start_radius-end_radius) < 1.e-6)
        {
            // Make an implicit truncated cylinder
            vtkSmartPointer<vtkCylinder> p_cyinder = vtkSmartPointer<vtkCylinder>::New();
            p_cyinder->SetRadius(start_radius);

            // Rotate the cylinder to align with the vessel vector
            c_vector<double, DIM> segment_vector = segments[idx]->GetUnitTangent();
            vtkSmartPointer<vtkTransform> p_transform = vtkSmartPointer<vtkTransform>::New();
            p_transform->PostMultiply();
            c_vector<double, DIM> mid_point = segments[idx]->GetMidPoint();

            p_transform->Translate(-mid_point[0], -mid_point[1], -mid_point[2]);
            double angle = std::acos(inner_prod(y_axis, segment_vector));
            if (std::abs(inner_prod(y_axis, segment_vector)) < 1.0 - 1.e-6)
            {
                c_vector<double, DIM> axis = VectorProduct(y_axis, segment_vector);
                p_transform->RotateWXYZ(-angle*180.0/M_PI, axis[0], axis[1], axis[2]);
            }
            p_cyinder->SetTransform(p_transform);

            // Create the clipping planes
            vtkSmartPointer<vtkPlane> p_start_plane = vtkSmartPointer<vtkPlane>::New();
            p_start_plane->SetOrigin(start_location[0], start_location[1], start_location[2]);
            p_start_plane->SetNormal(segment_vector[0], segment_vector[1], segment_vector[2]);

            vtkSmartPointer<vtkPlane> p_end_plane = vtkSmartPointer<vtkPlane>::New();
            p_end_plane->SetOrigin(end_location[0], end_location[1], end_location[2]);
            c_vector<double, DIM> opposite_segment_vector = -segment_vector;
            p_end_plane->SetNormal(opposite_segment_vector[0], opposite_segment_vector[1], opposite_segment_vector[2]);

            // Clip the cylinder
            vtkSmartPointer<vtkImplicitBoolean> p_truncated_cylinder = vtkSmartPointer<vtkImplicitBoolean>::New();
            p_truncated_cylinder->SetOperationTypeToDifference();
            p_truncated_cylinder->AddFunction(p_cyinder);
            p_truncated_cylinder->AddFunction(p_start_plane);
            p_truncated_cylinder->AddFunction(p_end_plane);

            // Add spheres to the end
            vtkSmartPointer<vtkSphere> p_start_sphere= vtkSmartPointer<vtkSphere>::New();
            p_start_sphere->SetRadius(start_radius);
            p_start_sphere->SetCenter(start_location[0], start_location[1], start_location[2]);

            vtkSmartPointer<vtkSphere> p_end_sphere= vtkSmartPointer<vtkSphere>::New();
            p_end_sphere->SetRadius(end_radius);
            p_end_sphere->SetCenter(end_location[0], end_location[1], end_location[2]);

            vtkSmartPointer<vtkImplicitBoolean> p_merged_segment = vtkSmartPointer<vtkImplicitBoolean>::New();
            p_merged_segment->SetOperationTypeToUnion();
            p_merged_segment->AddFunction(p_truncated_cylinder);
            p_merged_segment->AddFunction(p_start_sphere);
            p_merged_segment->AddFunction(p_end_sphere);

            segment_functions.push_back(p_merged_segment);
            //segment_functions.push_back(p_truncated_cylinder);
        }
        else if(start_radius > end_radius)
        {

        }
        else
        {

        }
    }

    // Reconstruct the surface
    vtkSmartPointer<vtkImplicitBoolean> p_merged_network = vtkSmartPointer<vtkImplicitBoolean>::New();
    p_merged_network->SetOperationTypeToUnion();

    if(segment_functions.size()>0)
    {
        p_merged_network->AddFunction(segment_functions[0]);
        if(segment_functions.size()>1)
        {
            for(unsigned idx=1; idx<segment_functions.size();idx++)
            {
                p_merged_network->AddFunction(segment_functions[idx]);
            }
        }
    }
//
    vtkSmartPointer<vtkSampleFunction> p_sampler = vtkSmartPointer<vtkSampleFunction>::New();
    p_sampler->SetImplicitFunction(p_merged_network);

    std::vector<std::pair<double, double> > extents = mpVesselNetwork->GetExtents();

    p_sampler->SetModelBounds(extents[0].first-100.0,
                              extents[0].second+100.0,
                              extents[1].first-100.0,
                              extents[1].second+100.0,
                              extents[2].first-100.0,
                              extents[2].second+100.0);
    p_sampler->SetSampleDimensions(100, 100, 100);
    vtkSmartPointer<vtkContourFilter> p_contour = vtkSmartPointer<vtkContourFilter>::New();
    p_contour->SetInputConnection(p_sampler->GetOutputPort());
    p_contour->SetValue(0, 0.0);
    p_contour->Update();

    mpSurface = p_contour->GetOutput();

    return mpSurface;
}

template class VesselNetworkMesher<3>;
