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

//#ifdef CHASTE_ANGIOGENESIS_VMTK
#include "Exception.hpp"
#include "VtkBoundaryExtractor.hpp"
#include <vtkFeatureEdges.h>
#include <vtkCleanPolyData.h>
#include <vtkTriangleFilter.h>
#include <vtkStripper.h>
#include <vtkSplineFilter.h>

VtkBoundaryExtractor::VtkBoundaryExtractor()
    : mpInputSurface(),
      mpOutputSurface(),
      mSmoothingLength(1),
      mDoSmoothing(true)
{

}

boost::shared_ptr<VtkBoundaryExtractor> VtkBoundaryExtractor::Create()
{
    MAKE_PTR(VtkBoundaryExtractor, pSelf);
    return pSelf;
}

VtkBoundaryExtractor::~VtkBoundaryExtractor()
{

}

vtkSmartPointer<vtkPolyData> VtkBoundaryExtractor::GetOutput()
{
    if(mpOutputSurface)
    {
        return(mpOutputSurface);
    }
    else
    {
        EXCEPTION("No output set. Did you run 'Update()' ?");
    }
}

void VtkBoundaryExtractor::SetInput(vtkSmartPointer<vtkPolyData> pInput)
{
    mpInputSurface = pInput;
}

void VtkBoundaryExtractor::SetSmoothingLength(double value)
{
    mSmoothingLength = value;
}

void VtkBoundaryExtractor::SetDoSmoothing(bool doSmoothing)
{
    mDoSmoothing = doSmoothing;
}

void VtkBoundaryExtractor::Update()
{
    if(!mpInputSurface)
    {
        EXCEPTION("No input set.");
    }

    vtkSmartPointer<vtkFeatureEdges> p_features = vtkSmartPointer<vtkFeatureEdges>::New();
    p_features->SetInput(mpInputSurface);
    p_features->Update();

    vtkSmartPointer<vtkCleanPolyData> p_clean = vtkSmartPointer<vtkCleanPolyData>::New();
    p_clean->SetInputConnection(p_features->GetOutputPort());
    p_clean->Update();

    vtkSmartPointer<vtkTriangleFilter> p_triangle = vtkSmartPointer<vtkTriangleFilter>::New();
    p_triangle->SetInputConnection(p_clean->GetOutputPort());
    p_triangle->Update();

    if(mDoSmoothing)
    {
        vtkSmartPointer<vtkStripper> p_stripper = vtkSmartPointer<vtkStripper>::New();
        p_stripper->SetInputConnection(p_triangle->GetOutputPort());
        p_stripper->Update();

        vtkSmartPointer<vtkSplineFilter> p_spline = vtkSmartPointer<vtkSplineFilter>::New();
        p_spline->SetInputConnection(p_stripper->GetOutputPort());
        p_spline->SetLength(mSmoothingLength);
        p_spline->Update();
        mpOutputSurface = p_spline->GetOutput();
    }
    else
    {
        mpOutputSurface = p_triangle->GetOutput();
    }
}
//#endif /*CHASTE_ANGIOGENESIS_VMTK*/
