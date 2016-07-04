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
#ifndef ImageIO_HPP_
#define ImageIO_HPP_

#include "SmartPointers.hpp"
#define _BACKWARD_BACKWARD_WARNING_H 1 //Cut out the vtk deprecated warning
#include <vtkImageData.h>
#include <vtkSmartPointer.h>
#include <itkImage.h>

/**
 *  This class manages input, output and conversion of 2D and 3D images, including VTK and ITK
 *  types.
 *
 *  Note: the supported image types are a small sub-set of those supported by VTK and ITK.
 */
class ImageIO
{

private:

    /**
     *  A vtk image
     */
    vtkSmartPointer<vtkImageData> mpVtkImage;

    /**
     *  An itk image. Not only unsigned char is used.
     */
    itk::Image<unsigned char, 2>::Pointer mpItkImage;

    /**
     *  A filepath
     */
    std::string mFilepath;

    /**
     *  The fraction of pixels to retain in X after reading
     */
    double mResizeX;

    /**
     *  The fraction of pixels to retain in Y after reading
     */
    double mResizeY;

    /**
     *  The fraction of pixels to retain in Z
     */
    double mResizeZ;

public:

    /**
     * Constructor
     */
    ImageIO();

    /**
     * Destructor
     */
    ~ImageIO();

    /**
     * Factory constructor method
     * @return a shared pointer to a instance of this class
     */
    static boost::shared_ptr<ImageIO> Create();

    /**
     * Convert data in mpVtkImage to mpItkImage. This involves changing the image to
     * unsigned char.
     */
    void ConvertVtkToItk();

    /**
     * Convert data in mpItkImage to mpVtkImage. This involves changing the image to
     * unsigned char.
     */
    void ConvertItkToVtk();

    /**
     * Get the image in vtk format
     */
    vtkSmartPointer<vtkImageData> GetVtkImage();

    /**
     * Get the image in itk format
     */
    itk::Image<unsigned char, 2>::Pointer GetItkImage();

    /**
     * Set the filename for the readers and writers
     * @param rFilename the file name
     */
    void SetFilename(const std::string& rFilename);

    /**
     * Set the resize factors for the image
     * @param factorX the fraction of pixels to retain in x
     * @param factorY the fraction of pixels to retain in y
     * @param factorZ the fraction of pixels to retain in z
     */
    void SetImageResizeFactors(double factorX, double factorY=1.0, double factorZ=1.0);

    /**
     * Set the image in vtk format
     */
    void SetVtkImage(vtkSmartPointer<vtkImageData> pImage);

    /**
     * Set the image in itk format
     */
    void SetItkImage(itk::Image<unsigned char, 2>::Pointer pImage);

    /**
     * Do the read operation
     */
    void ReadVtkImage();

    /**
     * Do the read operation
     */
    void ReadItkImage();

    /**
     * Write the image in VTK format
     */
    void WriteVtkImage();

    /**
     * Write the image in ITK format. TIFF or PNG only
     */
    void WriteItkImage();
};

#endif /*ImageIO_HPP_*/
#endif /*CHASTE_ANGIOGENESIS_EXTENDED*/
