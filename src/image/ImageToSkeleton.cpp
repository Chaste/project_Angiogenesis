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

#ifdef CHASTE_ANGIOGENESIS_EXTENDED
#include "Exception.hpp"
#include "ImageToSkeleton.hpp"
#include <itkImage.h>
#include <itkImageFileWriter.h>
#include <itkImageFileReader.h>
#include <itkBinaryThinningImageFilter.h>
#include <itkImageToVTKImageFilter.h>
#include <itkTIFFImageIOFactory.h>
#include <itkVTKImageToImageFilter.h>
#include <vtkImageSkeleton2D.h>
#include <vtkImageShiftScale.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkImageCast.h>
#include "itkPNGImageIOFactory.h"

ImageToSkeleton::ImageToSkeleton()
    : mpImage(vtkSmartPointer<vtkImageData>::New()),
      mpSkeleton(vtkSmartPointer<vtkImageData>::New()),
      mReverseIntensity(false),
      mUseVTKVersion(true)
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

void ImageToSkeleton::SetUseVtkVersion(bool value)
{
    mUseVTKVersion = value;
}

void ImageToSkeleton::Update()
{
    if(!mpImage)
    {
        EXCEPTION("No input set.");
    }

    if(mUseVTKVersion)
    {
        vtkSmartPointer<vtkImageShiftScale> imCast1 = vtkSmartPointer<vtkImageShiftScale>::New();
        imCast1->SetOutputScalarTypeToUnsignedChar();
        imCast1->SetInput(mpImage);
        imCast1->Update();

        vtkSmartPointer<vtkImageSkeleton2D> p_skeleton = vtkSmartPointer<vtkImageSkeleton2D>::New();
        p_skeleton->SetInput(imCast1->GetOutput());
        p_skeleton->SetNumberOfIterations(100);
        p_skeleton->SetPrune(0);
        p_skeleton->ReleaseDataFlagOff();
        p_skeleton->Update();
        mpSkeleton = p_skeleton->GetOutput();
    }
    else
    {
        itk::TIFFImageIOFactory::RegisterOneFactory();

        typedef itk::Image<unsigned char, 2> ImageType;
        typedef itk::ImageFileReader<ImageType> ImageReader;
        ImageReader::Pointer reader = ImageReader::New();
        reader->SetFileName("/home/grogan/median.tif");
        reader->Update();

//        vtkSmartPointer<vtkImageShiftScale> imCast1 = vtkSmartPointer<vtkImageShiftScale>::New();
//        imCast1->SetOutputScalarTypeToUnsignedChar();
//        imCast1->SetInput(mpImage);
//        imCast1->Update();

//        typedef itk::Image<unsigned char, 2> ImageType;
//        typedef itk::VTKImageToImageFilter<ImageType> VTKImageToImageType;
//        VTKImageToImageType::Pointer vtkImageToImageFilter = VTKImageToImageType::New();
//        vtkImageToImageFilter->SetInput(mpImage);
//        vtkImageToImageFilter->Update();

        typedef itk::BinaryThinningImageFilter<ImageType, ImageType> BinaryThinningImageFilterType;
        BinaryThinningImageFilterType::Pointer binaryThinningImageFilter = BinaryThinningImageFilterType::New();
//        binaryThinningImageFilter->SetInput(vtkImageToImageFilter->GetOutput());
        binaryThinningImageFilter->SetInput(reader->GetOutput());
        binaryThinningImageFilter->Update();

        typedef itk::ImageToVTKImageFilter<ImageType> ConnectorType;
        ConnectorType::Pointer connector = ConnectorType::New();
        connector->SetInput(binaryThinningImageFilter->GetOutput());
        connector->Update();
        connector->UpdateOutputInformation();

        vtkSmartPointer<vtkImageCast> imCast = vtkSmartPointer<vtkImageCast>::New();
        imCast->SetOutputScalarTypeToFloat();
        imCast->SetInput(connector->GetOutput());

        vtkImageData* p_output = imCast->GetOutput();
        p_output->Update();
        p_output->Register(NULL); // Avoids memory problems.
        mpSkeleton.TakeReference(p_output);

        // Todo move below to own class
//        itk::PNGImageIOFactory::RegisterOneFactory();
//        typedef  itk::ImageFileWriter< ImageType  > WriterType;
//        WriterType::Pointer writer = WriterType::New();
//        writer->SetFileName("/home/grogan/test.png");
//        writer->SetInput(vtkImageToImageFilter->GetOutput());
//        writer->Update();
    }

    // Convert to polydata. Iterate over the grid.



}
#endif /*CHASTE_ANGIOGENESIS_EXTENDED*/
