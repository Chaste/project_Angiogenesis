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
#ifndef ImageToSkeleton_HPP_
#define ImageToSkeleton_HPP_
#include "SmartPointers.hpp"
#include "VesselNetwork.hpp"

/**
* Get a vessel network from a binary mask in image format. The input image format should be an
* 8 bit grayscale TIFF. The (quick and dirty) 2D VTK algorithm or the considerably slower ITK binary thinning
* algorithms can be used. Output data in VTK formats at intermediate stages can be requested.
*/
class ImageToSkeleton
{
    /**
     *  Use the 2d vtk algorithm (default)
     */
    bool mUseVTK2d;

    /**
     *  Use the 2d itk algorithm
     */
    bool mUseITK2d;

    /**
     *  Use the 3d itk algorithm
     */
    bool mUseITK3d;

    /**
     *  Verbose output
     */
    bool mVerboseOutput;

    /**
     *  Prune level for vtk algorithm
     */
    unsigned mVtkPruneLevel;

    /**
     *  The output directory
     */
    std::string mOutputDirectory;

    /**
     *  The input file
     */
    std::string mInputFile;

    /**
     *  The vessel network
     */
    boost::shared_ptr<VesselNetwork<3> > mpVesselNetwork;

public:

    /**
     *  Constructor
     */
    ImageToSkeleton();

    /**
     *  Destructor
     */
    ~ImageToSkeleton();

    /**
     *  Factory constructor method
     */
    static boost::shared_ptr<ImageToSkeleton> Create();

    /**
     *  Return a vessel network
     */
    boost::shared_ptr<VesselNetwork<3> > GetOutput();

    /**
     *  Set filename
     */
    void SetFilename(std::string filename);

    /**
     *  Set the output directory for verbose mode
     */
    void SetOutputDirectory(std::string directory);

    /**
     *  Use the vtk skeletonize algorithm
     */
    void SetUseVtkAlgorithm(bool value);

    /**
     *  Use the itk 2d binary thinning algorithm
     */
    void SetUseItk2dAlgorithm(bool value);

    /**
     *  Use the itk 3d binary thinning algorithm
     */
    void SetUseItk3dAlgorithm(bool value);

    /**
     *  Do verbose output
     */
    void SetUseVerboseMode(bool value);

    /**
     *  Run the tool
     */
    void Update();

};

#endif /*ImageToSkeleton_HPP_*/
#endif /*CHASTE_ANGIOGENESIS_EXTENDED*/
