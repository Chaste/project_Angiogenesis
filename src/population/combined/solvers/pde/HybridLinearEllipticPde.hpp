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

#ifndef HYBRIDLINEARELLIPTICPDE_HPP_
#define HYBRIDLINEARELLIPTICPDE_HPP_

#include <string>
#define _BACKWARD_BACKWARD_WARNING_H 1 //Cut out the vtk deprecated warning for now (gcc4.3)
#include <vtkBox.h>
#include <vtkTetra.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCellLocator.h>
#include "ChastePoint.hpp"
#include "UblasIncludes.hpp"
#include "SmartPointers.hpp"
#include "UblasVectorInclude.hpp"
#include "AbstractLinearEllipticPde.hpp"
#include "SimpleCellPopulation.hpp"
#include "CaVascularNetwork.hpp"

/*
 * Linear Elliptic PDE with discrete or averaged cell and
 * vessel source terms.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM = ELEMENT_DIM>
class HybridLinearEllipticPde : public AbstractLinearEllipticPde<ELEMENT_DIM, SPACE_DIM>
{
    c_matrix<double, SPACE_DIM, SPACE_DIM> mDiffusionTensor;
    double mDiffusivity;
    double mConstantInUTerm;
    double mLinearInUTerm;
    std::string mVariableName;
    boost::shared_ptr<SimpleCellPopulation> mpPopulation;
    boost::shared_ptr<CaVascularNetwork<3> > mpNetwork;

public:

    HybridLinearEllipticPde() :
            mDiffusionTensor(identity_matrix<double>(SPACE_DIM)),
            mDiffusivity(1.e-3),
            mConstantInUTerm(0.0),
            mLinearInUTerm(0.0),
            mVariableName("Default"),
            mpPopulation(),
            mpNetwork()
    {
        mDiffusionTensor *= mDiffusivity;
    }

    static boost::shared_ptr<HybridLinearEllipticPde<ELEMENT_DIM, SPACE_DIM> > Create()
    {
        MAKE_PTR(HybridLinearEllipticPde<ELEMENT_DIM>, pSelf);
        return pSelf;
    }

    double ComputeConstantInUSourceTerm(const ChastePoint<SPACE_DIM>& rX, Element<ELEMENT_DIM, SPACE_DIM>* pElement)
    {
        return mConstantInUTerm;
    }

    double ComputeLinearInUCoeffInSourceTerm(const ChastePoint<SPACE_DIM>& rX, Element<ELEMENT_DIM, SPACE_DIM>* pElement)
    {
        return mLinearInUTerm;
    }

    c_matrix<double, SPACE_DIM, SPACE_DIM> ComputeDiffusionTerm(const ChastePoint<SPACE_DIM>&)
    {
        return mDiffusionTensor;
    }

    void SetCellPopulation(boost::shared_ptr<SimpleCellPopulation> pPopulation)
    {
        mpPopulation = pPopulation;
    }

    void SetVesselNetwork(boost::shared_ptr<CaVascularNetwork<3> > pNetwork)
    {
        mpNetwork = pNetwork;
    }

    void SetConstantInUTerm(double constantInUTerm)
    {
        mConstantInUTerm = constantInUTerm;
    }

    void SetLinearInUTerm(double linearInUTerm)
    {
        mLinearInUTerm = linearInUTerm;
    }

    void SetDiffusionConstant(double diffusivity)
    {
        mDiffusivity = diffusivity;
        mDiffusionTensor = identity_matrix<double>(SPACE_DIM)* mDiffusivity;
    }

    double GetConstantInUTerm(c_vector<double, 3> location = zero_vector<double>(3), double spacing = 0.0)
    {
        return mConstantInUTerm;
    }

    void SetVariableName(const std::string& rVariableName)
    {
        mVariableName = rVariableName;
    }

    const std::string& GetVariableName()
    {
        return mVariableName;
    }

    double GetLinearInUTerm(c_vector<double, 3> location  = zero_vector<double>(3), double spacing = 0.0)
    {
        double cell_consumption_term = 0.0;
        if(mpPopulation)
        {
            unsigned num_points = 0;
            for (unsigned mdx = 0; mdx < mpPopulation->GetCells().size(); mdx++)
            {
                c_vector<double, 3> cell_location = mpPopulation->GetCells()[mdx]->rGetLocation();
                if (IsPointInBox(cell_location, location, spacing))
                {
                    num_points++;
                }
            }
            cell_consumption_term = 1.e-5 * double(num_points);
        }
        return mLinearInUTerm - cell_consumption_term;
    }

    double GetDiffusionConstant()
    {
        return mDiffusivity;
    }

    bool IsPointInBox(c_vector<double, 3> point, c_vector<double, 3> location, double spacing)
    {
        bool point_in_box = false;
        if(point[0] >= location[0] -spacing/2.0 && point[0] <= location [0] + spacing/2.0)
        {
            if(point[1] >= location[1] -spacing/2.0 && point[1] <= location [1] + spacing/2.0)
            {
                if(point[2] >= location[2] -spacing/2.0 && point[2] <= location [2] + spacing/2.0)
                {
                    return true;
                }
            }
        }
        return point_in_box;
    }

    bool IsPointInTetra(c_vector<double, 3> start_point, std::vector<c_vector<double, 3> > locations)
    {
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

    double LengthOfLineInBox(c_vector<double, 3> start_point, c_vector<double, 3> end_point, c_vector<double, 3> location, double spacing)
    {
        // If the line is fully in the box return its length
        bool point1_in_box = IsPointInBox(start_point, location, spacing);
        bool point2_in_box = IsPointInBox(end_point, location, spacing);
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
            bounds[4] = location[2] - spacing/2.0;
            bounds[5] = location[2] + spacing/2.0;

            double t1;
            double t2;
            int plane1;
            int plane2;
            c_vector<double,3> intercept_1;
            c_vector<double,3> intercept_2;

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

    double LengthOfLineInTetra(c_vector<double, 3> start_point, c_vector<double, 3> end_point, std::vector<c_vector<double, 3> > locations)
    {
        bool point1_in_tetra = IsPointInTetra(start_point, locations);
        bool point2_in_tetra = IsPointInTetra(end_point, locations);

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
            c_vector<double,3> intersection;
            c_vector<double,3> parametric_intersection;
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
                c_vector<double,3> intersection2;
                p_grid->GetCell(0)->IntersectWithLine(&end_point[0], &start_point[0], 1.e-6, t, &intersection2[0], &parametric_intersection[0], subId);
                return norm_2(intersection - intersection2);
            }
            else
            {
                return 0.0;
            }
        }
    }
};

#endif /*HYBRIDLINEARELLIPTICPDE_HPP_*/
