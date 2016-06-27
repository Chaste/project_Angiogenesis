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
#include "ImageToSkeleton.hpp"
#include <vtkImageSkeleton2D.h>
#include <vtkPointData.h>
#include <vtkCellData.h>

ImageToSkeleton::ImageToSkeleton()
    : mpImage(),
      mpSkeleton(),
      mReverseIntensity(false)
{

}

boost::shared_ptr<ImageToSkeleton> ImageToSkeleton::Create()
{
    MAKE_PTR(ImageToSkeleton, pSelf);
    return pSelf;
}

ImageToSkeleton::~ImageToSkeleton()
{

}

void ImageToSkeleton::SetReverseIntensity(bool value)
{
    mReverseIntensity = value;
}

vtkSmartPointer<vtkImageData> ImageToSkeleton::GetOutput()
{
    if(mpSkeleton)
    {
        return(mpSkeleton);
    }
    else
    {
        EXCEPTION("No output set. Did you run 'Update()' ?");
    }
}

void ImageToSkeleton::SetInput(vtkSmartPointer<vtkImageData> pImage)
{
    mpImage = pImage;
}

void ImageToSkeleton::Update()
{
    if(!mpImage)
    {
        EXCEPTION("No input set.");
    }

    if(mReverseIntensity)
    {
        double max_intensity = mpImage->GetPointData()->GetArray("ImageScalars")->GetDataTypeMax();
        if(max_intensity>0.0)
        {
            for(unsigned idx=0; idx< mpImage->GetNumberOfPoints(); idx++)
            {
                double current_intensity = mpImage->GetPointData()->GetArray("ImageScalars")->GetTuple1(idx);
                double new_value = max_intensity - current_intensity;
                if(new_value==1.0)
                {
                    new_value = 2.0;
                }
                mpImage->GetPointData()->GetArray("ImageScalars")->SetTuple1(idx, new_value);
            }
        }
    }
    mpImage->GetPointData()->SetScalars(mpImage->GetPointData()->GetArray("ImageScalars"));
    mpImage->GetCellData()->SetScalars(mpImage->GetPointData()->GetArray("ImageScalars"));

    vtkSmartPointer<vtkImageSkeleton2D> p_skeleton = vtkSmartPointer<vtkImageSkeleton2D>::New();
    p_skeleton->SetInput(mpImage);
    p_skeleton->SetNumberOfIterations(20);
    p_skeleton->ReleaseDataFlagOff();
    p_skeleton->Update();
    mpSkeleton = p_skeleton->GetOutput();

}
//#endif /*CHASTE_ANGIOGENESIS_VMTK*/
