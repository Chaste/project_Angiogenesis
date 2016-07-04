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
#include <boost/filesystem.hpp>
#include <vtkTIFFReader.h>
#include <vtkPNGReader.h>
#include <vtkJPEGReader.h>
#include <vtkBMPReader.h>
#include <vtkImageResize.h>
#include <vtkImageLuminance.h>
#include <vtkImageCast.h>
#include <vtkXMLImageDataWriter.h>
#include <itkVTKImageToImageFilter.h>
#include <itkTIFFImageIOFactory.h>
#include <itkPNGImageIOFactory.h>
#include <itkImageFileReader.h>
#include <itkImageToVTKImageFilter.h>
#include <itkImageFileWriter.h>
#include "ImageIO.hpp"

ImageIO::ImageIO()
    : mpVtkImage(vtkSmartPointer<vtkImageData>::New()),
      mpItkImage(itk::Image<unsigned char, 2>::New()),
      mFilepath(),
      mResizeX(1.0),
      mResizeY(1.0),
      mResizeZ(1.0)
{

}

boost::shared_ptr<ImageIO> ImageIO::Create()
{
    MAKE_PTR(ImageIO, pSelf);
    return pSelf;
}

ImageIO::~ImageIO()
{

}

void ImageIO::ConvertVtkToItk()
{
    if(!mpVtkImage)
    {
        EXCEPTION("Vtk image data has not been set.");
    }

    vtkSmartPointer<vtkImageCast> imCast = vtkSmartPointer<vtkImageCast>::New();
    imCast->SetOutputScalarTypeToUnsignedChar();
    imCast->SetInput(mpVtkImage);
    imCast->Update();

    typedef itk::Image<unsigned char, 2> ImageType;
    typedef itk::VTKImageToImageFilter<ImageType> VTKImageToImageType;
    VTKImageToImageType::Pointer vtkImageToImageFilter = VTKImageToImageType::New();
    vtkImageToImageFilter->SetInput(imCast->GetOutput());
    vtkImageToImageFilter->Update();
    mpItkImage = vtkImageToImageFilter->GetOutput();
}

void ImageIO::ConvertItkToVtk()
{
    if(!mpItkImage)
    {
        EXCEPTION("Itk image data has not been set.");
    }

    typedef itk::Image<unsigned char, 2> ImageType;
    typedef itk::ImageToVTKImageFilter<ImageType> ConnectorType;
    ConnectorType::Pointer connector = ConnectorType::New();
    connector->SetInput(mpItkImage);
    connector->Update();
    connector->UpdateOutputInformation();

    vtkSmartPointer<vtkImageCast> imCast = vtkSmartPointer<vtkImageCast>::New();
    imCast->SetOutputScalarTypeToFloat();
    imCast->SetInput(connector->GetOutput());

    vtkImageData* p_output = imCast->GetOutput();
    p_output->Update();
    p_output->Register(NULL); // Avoids memory problems.
    mpVtkImage.TakeReference(p_output);
}

vtkSmartPointer<vtkImageData> ImageIO::GetVtkImage()
{
    if(mpVtkImage)
    {
        return mpVtkImage;
    }
    else
    {
        if(mpItkImage)
        {
            // Convert to VTK
            ConvertItkToVtk();
            return mpVtkImage;
        }
        else
        {
            EXCEPTION("No image data has been set.");
        }
    }
}

itk::Image<unsigned char, 2>::Pointer ImageIO::GetItkImage()
{
    if(mpItkImage)
    {
        return mpItkImage;
    }
    else
    {
        if(mpVtkImage)
        {
            // Convert to ITK
            ConvertVtkToItk();
            return mpItkImage;
        }
        else
        {
            EXCEPTION("No image data has been set.");
        }
    }
}

void ImageIO::SetFilename(const std::string& filename)
{
    mFilepath = filename;
}

void ImageIO::SetImageResizeFactors(double factorX, double factorY, double factorZ)
{
    mResizeX = factorX;
    mResizeY = factorY;
    mResizeZ = factorZ;
}

void ImageIO::ReadVtkImage()
{
    // Get the file extension
    if(mFilepath == "")
    {
        EXCEPTION("Input file not specified for image reader");
    }
    std::string file_extension  = boost::filesystem::extension(mFilepath);

    if(file_extension == ".tif" or file_extension == ".tiff" or file_extension == ".TIF" or file_extension == ".TIFF")
    {
        vtkSmartPointer<vtkTIFFReader> p_reader = vtkSmartPointer<vtkTIFFReader>::New();
        p_reader->SetFileName(mFilepath.c_str());
        if(mResizeX==1.0 and mResizeY==1.0 and mResizeZ==1.0)
        {
            p_reader->Update();
            mpVtkImage = p_reader->GetOutput();
        }
        else
        {
            vtkSmartPointer<vtkImageResize> p_resize = vtkSmartPointer<vtkImageResize>::New();
            p_resize->SetInputConnection(p_reader->GetOutputPort());
            p_resize->SetMagnificationFactors(mResizeX, mResizeY, mResizeZ);
            p_resize->SetResizeMethodToMagnificationFactors();
            p_resize->Update();
            mpVtkImage = p_resize->GetOutput();
        }
    }
    else if(file_extension == ".png" or file_extension == ".PNG")
    {
        vtkSmartPointer<vtkPNGReader> p_reader = vtkSmartPointer<vtkPNGReader>::New();
        p_reader->SetFileName(mFilepath.c_str());
        if(mResizeX==1.0 and mResizeY==1.0 and mResizeZ==1.0)
        {
            p_reader->Update();

            // RGB to scalar
            vtkSmartPointer<vtkImageLuminance> p_lum = vtkSmartPointer<vtkImageLuminance>::New();
            p_lum->SetInputConnection(p_reader->GetOutputPort());
            p_lum->Update();
            mpVtkImage = p_lum->GetOutput();
        }
        else
        {
            vtkSmartPointer<vtkImageResize> p_resize = vtkSmartPointer<vtkImageResize>::New();
            p_resize->SetInputConnection(p_reader->GetOutputPort());
            p_resize->SetMagnificationFactors(mResizeX, mResizeY, mResizeZ);
            p_resize->SetResizeMethodToMagnificationFactors();
            p_resize->Update();

            // RGB to scalar
            vtkSmartPointer<vtkImageLuminance> p_lum = vtkSmartPointer<vtkImageLuminance>::New();
            p_lum->SetInputConnection(p_resize->GetOutputPort());
            p_lum->Update();
            mpVtkImage = p_lum->GetOutput();
        }
    }
    else if(file_extension == ".jpg" or file_extension == ".JPG")
    {
        vtkSmartPointer<vtkJPEGReader> p_reader = vtkSmartPointer<vtkJPEGReader>::New();
        p_reader->SetFileName(mFilepath.c_str());
        if(mResizeX==1.0 and mResizeY==1.0 and mResizeZ==1.0)
        {
            p_reader->Update();
            mpVtkImage = p_reader->GetOutput();
        }
        else
        {
            vtkSmartPointer<vtkImageResize> p_resize = vtkSmartPointer<vtkImageResize>::New();
            p_resize->SetInputConnection(p_reader->GetOutputPort());
            p_resize->SetMagnificationFactors(mResizeX, mResizeY, mResizeZ);
            p_resize->SetResizeMethodToMagnificationFactors();
            p_resize->Update();
            mpVtkImage = p_resize->GetOutput();
        }
    }
    else if(file_extension == ".bmp" or file_extension == ".BMP")
    {
        vtkSmartPointer<vtkBMPReader> p_reader = vtkSmartPointer<vtkBMPReader>::New();
        p_reader->SetFileName(mFilepath.c_str());
        if(mResizeX==1.0 and mResizeY==1.0 and mResizeZ==1.0)
        {
            p_reader->Update();
            mpVtkImage = p_reader->GetOutput();
        }
        else
        {
            vtkSmartPointer<vtkImageResize> p_resize = vtkSmartPointer<vtkImageResize>::New();
            p_resize->SetInputConnection(p_reader->GetOutputPort());
            p_resize->SetMagnificationFactors(mResizeX, mResizeY, mResizeZ);
            p_resize->SetResizeMethodToMagnificationFactors();
            p_resize->Update();
            mpVtkImage = p_resize->GetOutput();
        }
    }
    else
    {
        EXCEPTION("Input file extension not recognized");
    }

    if(!mpVtkImage)
    {
        EXCEPTION("Image reading failed.");
    }
}

void ImageIO::ReadItkImage()
{
    // Get the file extension
    if(mFilepath == "")
    {
        EXCEPTION("Input file not specified for image reader");
    }
    std::string file_extension  = boost::filesystem::extension(mFilepath);

    if(file_extension == ".tif" or file_extension == ".tiff" or file_extension == ".TIF" or file_extension == ".TIFF")
    {
        itk::TIFFImageIOFactory::RegisterOneFactory();
        typedef itk::Image<unsigned char, 2> ImageType;
        typedef itk::ImageFileReader<ImageType> ImageReader;
        ImageReader::Pointer reader = ImageReader::New();
        reader->SetFileName(mFilepath.c_str());
        reader->Update();
        mpItkImage = reader->GetOutput();
    }
    else if(file_extension == ".png" or file_extension == ".PNG")
    {
        itk::PNGImageIOFactory::RegisterOneFactory();
        typedef itk::Image<unsigned char, 2> ImageType;
        typedef itk::ImageFileReader<ImageType> ImageReader;
        ImageReader::Pointer reader = ImageReader::New();
        reader->SetFileName(mFilepath.c_str());
        reader->Update();
        mpItkImage = reader->GetOutput();
    }
    else
    {
        EXCEPTION("Input file extension not supported");
    }

    if(!mpItkImage)
    {
        EXCEPTION("Image reading failed.");
    }
}

void ImageIO::WriteVtkImage()
{
    if(mFilepath == "")
    {
        EXCEPTION("Output file not specified for image writer");
    }

    vtkSmartPointer<vtkXMLImageDataWriter> p_writer1 = vtkSmartPointer<vtkXMLImageDataWriter>::New();
    p_writer1->SetFileName(mFilepath.c_str());
    p_writer1->SetInput(GetVtkImage());
    p_writer1->Write();
}

void ImageIO::WriteItkImage()
{
    if(mFilepath == "")
    {
        EXCEPTION("Output file not specified for image writer");
    }

    itk::TIFFImageIOFactory::RegisterOneFactory();
    itk::PNGImageIOFactory::RegisterOneFactory();
    typedef itk::Image<unsigned char, 2> ImageType;
    typedef itk::ImageFileWriter<ImageType> WriterType;
    WriterType::Pointer writer = WriterType::New();
    writer->SetFileName(mFilepath.c_str());
    writer->SetInput(GetItkImage());
    writer->Update();
}

#endif /*CHASTE_ANGIOGENESIS_EXTENDED*/
