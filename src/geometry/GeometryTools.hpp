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

#ifndef GEOMETRYTOOLS_HPP_
#define GEOMETRYTOOLS_HPP_

#include <vector>
#define _BACKWARD_BACKWARD_WARNING_H 1 //Cut out the vtk deprecated warning for now (gcc4.3)
#include <vtkBox.h>
#include <vtkTetra.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCellLocator.h>
#include "ChastePoint.hpp"
#include "SmartPointers.hpp"
#include "UblasVectorInclude.hpp"

/*
 * A collection of useful geometry functions
 */

/*
 * Is the point inside the cube box defined by a centre location and box side length
 */
template<unsigned DIM>
bool IsPointInBox(c_vector<double, DIM> point, c_vector<double, DIM> location, double spacing)
{
    bool point_in_box = false;
    if(point[0] >= location[0] -spacing/2.0 && point[0] <= location [0] + spacing/2.0)
    {
        if(point[1] >= location[1] -spacing/2.0 && point[1] <= location [1] + spacing/2.0)
        {
            if(DIM==3)
            {
                if(point[2] >= location[2] -spacing/2.0 && point[2] <= location [2] + spacing/2.0)
                {
                    return true;
                }
            }
            else
            {
                return true;
            }
        }
    }
    return point_in_box;
}

/*
 * Is the point inside the tetrahedron given by the vector of vertex locations
 */
template<unsigned DIM>
bool IsPointInTetra(c_vector<double, DIM> start_point, std::vector<c_vector<double, DIM> > locations)
{
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints> :: New();
    if(DIM==3)
    {
        points->InsertNextPoint(locations[0][0], locations[0][1], locations[0][2]);
        points->InsertNextPoint(locations[1][0], locations[1][1], locations[1][2]);
        points->InsertNextPoint(locations[2][0], locations[2][1], locations[2][2]);
        points->InsertNextPoint(locations[3][0], locations[3][1], locations[3][2]);
    }
    else
    {
        points->InsertNextPoint(locations[0][0], locations[0][1], 0.0);
        points->InsertNextPoint(locations[1][0], locations[1][1], 0.0);
        points->InsertNextPoint(locations[2][0], locations[2][1], 0.0);
        points->InsertNextPoint(locations[3][0], locations[3][1], 0.0);
    }


    vtkSmartPointer<vtkUnstructuredGrid> p_grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    p_grid->SetPoints(points);

    vtkIdType ptIds[] = {0, 1, 2, 3};
    p_grid->InsertNextCell(VTK_TETRA, 4, ptIds);
    p_grid->Update();

    vtkSmartPointer<vtkCellLocator> p_locator = vtkSmartPointer<vtkCellLocator>::New();
    p_locator->SetDataSet(p_grid);
    p_locator->Update();
    int in_tetra = p_locator->FindCell(&start_point[0]);
    if(in_tetra == -1)
    {
        return false;
    }
    else
    {
        return true;
    }
}

/*
 * Return the length of the line given by a start point and end point in the box given by a centre location and side length
 */
template<unsigned DIM>
double LengthOfLineInBox(c_vector<double, DIM> start_point,
                         c_vector<double, DIM> end_point,
                         c_vector<double, DIM> location, double spacing)
{
    // If the line is fully in the box return its length
    bool point1_in_box = IsPointInBox<DIM>(start_point, location, spacing);
    bool point2_in_box = IsPointInBox<DIM>(end_point, location, spacing);
    if(point1_in_box && point2_in_box)
    {
        return norm_2(end_point - start_point);
    }
    else
    {
        c_vector<double,6> bounds;
        bounds[0] = location[0] - spacing/2.0;
        bounds[1] = location[0] + spacing/2.0;
        bounds[2] = location[1] - spacing/2.0;
        bounds[3] = location[1] + spacing/2.0;
        if(DIM==3)
        {
            bounds[4] = location[2] - spacing/2.0;
            bounds[5] = location[2] + spacing/2.0;
        }
        else
        {
            bounds[4] = 0.0;
            bounds[5] = 0.0;
        }

        double t1;
        double t2;
        int plane1;
        int plane2;
        c_vector<double,DIM> intercept_1;
        c_vector<double,DIM> intercept_2;

        int in_box = vtkBox::IntersectWithLine(&bounds[0], &start_point[0], &end_point[0], t1, t2, &intercept_1[0], &intercept_2[0], plane1, plane2);

        if(point1_in_box)
        {
            return norm_2(intercept_2 - start_point);
        }

        if(point2_in_box)
        {
            return norm_2(intercept_1 - end_point);
        }

        if(in_box)
        {
            return norm_2(intercept_2 - intercept_1);
        }
        else
        {
            return 0.0;
        }
    }
}

/*
 * Return the length of the line given by a start point and end point in the tetrahedron given by vertex locations
 */
template<unsigned DIM>
double LengthOfLineInTetra(c_vector<double, DIM> start_point,
                           c_vector<double, DIM> end_point,
                           std::vector<c_vector<double, DIM> > locations)
{
    bool point1_in_tetra = IsPointInTetra<DIM>(start_point, locations);
    bool point2_in_tetra = IsPointInTetra<DIM>(end_point, locations);

    if(point1_in_tetra && point2_in_tetra)
    {
        return norm_2(end_point - start_point);
    }
    else
    {
        int line_crosses;

        vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints> :: New();
        points->InsertNextPoint(locations[0][0], locations[0][1], locations[0][2]);
        points->InsertNextPoint(locations[1][0], locations[1][1], locations[1][2]);
        points->InsertNextPoint(locations[2][0], locations[2][1], locations[2][2]);
        points->InsertNextPoint(locations[3][0], locations[3][1], locations[3][2]);

        vtkSmartPointer<vtkUnstructuredGrid> p_grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
        p_grid->SetPoints(points);

        vtkIdType ptIds[] = {0, 1, 2, 3};
        p_grid->InsertNextCell(VTK_TETRA, 4, ptIds);
        p_grid->Update();

        double t;
        c_vector<double,DIM> intersection;
        c_vector<double,DIM> parametric_intersection;
        int subId;

        if(point1_in_tetra)
        {
            p_grid->GetCell(0)->IntersectWithLine(&start_point[0], &end_point[0], 1.e-6, t, &intersection[0], &parametric_intersection[0], subId);
            return norm_2(intersection - start_point);
        }

        if(point2_in_tetra)
        {
            p_grid->GetCell(0)->IntersectWithLine(&end_point[0], &start_point[0], 1.e-6, t, &intersection[0], &parametric_intersection[0], subId);
            return norm_2(intersection - end_point);
        }

        line_crosses = p_grid->GetCell(0)->IntersectWithLine(&start_point[0], &end_point[0], 1.e-6, t, &intersection[0], &parametric_intersection[0], subId);
        if(line_crosses)
        {
            c_vector<double,DIM> intersection2;
            p_grid->GetCell(0)->IntersectWithLine(&end_point[0], &start_point[0], 1.e-6, t, &intersection2[0], &parametric_intersection[0], subId);
            return norm_2(intersection - intersection2);
        }
        else
        {
            return 0.0;
        }
    }
}

/*
 * Rotate the supplied vector about the axis by the specified angle.
 */
template<unsigned DIM>
c_vector<double, DIM> RotateAboutAxis(c_vector<double, DIM> direction, c_vector<double, DIM> axis, double angle)
{
    double sin_a = std::sin(angle);
    double cos_a = std::cos(angle);
    c_vector<double, DIM> unit_axis = axis / norm_2(axis);

    double dot_product = inner_prod(direction, unit_axis);
    c_vector<double, DIM> new_direction;

    if(DIM==3)
    {
        new_direction[0] = (unit_axis[0] * dot_product * (1.0 - cos_a) + direction[0] * cos_a
                    + (-unit_axis[2] * direction[1] + unit_axis[1] * direction[2]) * sin_a);
        new_direction[1] = (unit_axis[1] * dot_product * (1.0 - cos_a) + direction[1] * cos_a
                    + (unit_axis[2] * direction[0] - unit_axis[0] * direction[2]) * sin_a);
        new_direction[2] = (unit_axis[2] * dot_product * (1.0 - cos_a) + direction[2] * cos_a
                    + (-unit_axis[1] * direction[0] + unit_axis[0] * direction[1]) * sin_a);
    }
    else
    {
        new_direction[0] = unit_axis[0] * dot_product * (1.0 - cos_a) + direction[0] * cos_a;
        new_direction[1] = unit_axis[1] * dot_product * (1.0 - cos_a) + direction[1] * cos_a;
    }
    return new_direction;
}

#endif /*GEOMETRYTOOLS_HPP_*/
