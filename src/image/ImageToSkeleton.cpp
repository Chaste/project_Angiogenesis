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
#include <math.h>
#define _BACKWARD_BACKWARD_WARNING_H 1 //Cut out the vtk deprecated warning
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkBinaryThinningImageFilter.h>
//#include <itkBinaryThinningImageFilter3D.h>
#include <itkImageToVTKImageFilter.h>
#include <itkTIFFImageIOFactory.h>
#include <vtkImageData.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkImageSkeleton2D.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkContourFilter.h>
#include <vtkKdTreePointLocator.h>
#include <vtkImageCast.h>
#include "ImageIO.hpp"
#include "ImageToSkeleton.hpp"

ImageToSkeleton::ImageToSkeleton()
    : mUseVTK2d(true),
      mUseITK2d(false),
      mUseITK3d(false),
      mVerboseOutput(true),
      mVtkPruneLevel(0),
      mOutputDirectory(),
      mInputFile(),
      mpVesselNetwork()

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

boost::shared_ptr<VesselNetwork<3> > ImageToSkeleton::GetOutput()
{
    if(mpVesselNetwork)
    {
        return(mpVesselNetwork);
    }
    else
    {
        EXCEPTION("No output set. Did you run 'Update()' ?");
    }
}


void ImageToSkeleton::SetFilename(std::string filename)
{
    mInputFile = filename;
}

void ImageToSkeleton::SetOutputDirectory(std::string directory)
{
    mOutputDirectory = directory;
}

void ImageToSkeleton::SetUseVtkAlgorithm(bool value)
{
    mUseVTK2d = value;
}

void ImageToSkeleton::SetUseItk2dAlgorithm(bool value)
{
    mUseITK2d = value;
}

void ImageToSkeleton::SetUseItk3dAlgorithm(bool value)
{
    mUseITK3d = value;
}

void ImageToSkeleton::Update()
{
//    if(mInputFile.empty())
//    {
//        EXCEPTION("No input file set.");
//    }
//
//    vtkSmartPointer<vtkImageData> p_image;
//    if(mUseVTK2d)
//    {
//        // Read the image
//        ImageIO reader;
//        reader.SetFilename(mInputFile);
//        reader.Read()
//
//        vtkSmartPointer<vtkImageSkeleton2D> p_skeleton = vtkSmartPointer<vtkImageSkeleton2D>::New();
//        p_skeleton->SetInput(reader.GetOutput());
//        p_skeleton->SetNumberOfIterations(300);
//        p_skeleton->SetPrune(mVtkPruneLevel);
//        p_skeleton->ReleaseDataFlagOff();
//        p_skeleton->Update();
//
//        p_image = p_skeleton->GetOutput();
//    }
//    else if(mUseITK2d)
//    {
//        itk::TIFFImageIOFactory::RegisterOneFactory();
//        typedef itk::Image<unsigned char, 2> ImageType;
//        typedef itk::ImageFileReader<ImageType> ImageReader;
//        ImageReader::Pointer reader = ImageReader::New();
//        reader->SetFileName(mInputFile);
//        reader->Update();
//
//        typedef itk::BinaryThinningImageFilter<ImageType, ImageType> BinaryThinningImageFilterType;
//        BinaryThinningImageFilterType::Pointer binaryThinningImageFilter = BinaryThinningImageFilterType::New();
//        binaryThinningImageFilter->SetInput(reader->GetOutput());
//        binaryThinningImageFilter->Update();
//
//        typedef itk::ImageToVTKImageFilter<ImageType> ConnectorType;
//        ConnectorType::Pointer connector = ConnectorType::New();
//        connector->SetInput(binaryThinningImageFilter->GetOutput());
//        connector->Update();
//        connector->UpdateOutputInformation();
//
//        vtkSmartPointer<vtkImageCast> imCast = vtkSmartPointer<vtkImageCast>::New();
//        imCast->SetOutputScalarTypeToFloat();
//        imCast->SetInput(connector->GetOutput());
//
//        vtkImageData* p_output = imCast->GetOutput();
//        p_output->Update();
//        p_output->Register(NULL); // Avoids memory problems.
//        p_image.TakeReference(p_output);
//    }
//    else if(mUseITK3d)
//    {
//        itk::TIFFImageIOFactory::RegisterOneFactory();
//        typedef itk::Image<unsigned char, 3> ImageType;
//        typedef itk::ImageFileReader<ImageType> ImageReader;
//        ImageReader::Pointer reader = ImageReader::New();
//        reader->SetFileName(mInputFile);
//        reader->Update();
//
//        typedef itk::BinaryThinningImageFilter3D<ImageType, ImageType> BinaryThinningImageFilter3DType;
//        BinaryThinningImageFilter3DType::Pointer binaryThinningImageFilter3d = BinaryThinningImageFilter3DType::New();
//        binaryThinningImageFilter3d->SetInput(reader->GetOutput());
//        binaryThinningImageFilter3d->Update();
//
//        typedef itk::ImageToVTKImageFilter<ImageType> ConnectorType;
//        ConnectorType::Pointer connector = ConnectorType::New();
//        connector->SetInput(binaryThinningImageFilter->GetOutput());
//        connector->Update();
//        connector->UpdateOutputInformation();
//
//        vtkSmartPointer<vtkImageCast> imCast = vtkSmartPointer<vtkImageCast>::New();
//        imCast->SetOutputScalarTypeToFloat();
//        imCast->SetInput(connector->GetOutput());
//
//        vtkImageData* p_output = imCast->GetOutput();
//        p_output->Update();
//        p_output->Register(NULL); // Avoids memory problems.
//        p_image.TakeReference(p_output);
//    }
//    else
//    {
//        EXCEPTION("Skeleton algorithm not chosen")
//    }
//
//    if(!p_image)
//    {
//        EXCEPTION("Skeleton extraction failed in external skeletonize algorithms.");
//    }
//
//    if(mVerboseOutput and !mOutputDirectory.empty())
//    {
//        ImageIO writer;
//        writer.SetFilename(mOutputDirectory+"/skeleton.vti");
//        writer.SetImage(p_image);
//        writer.Write();
//    }
//
//    // Convert the pixelated skeleton into polylines
//    vtkSmartPointer<vtkContourFilter> p_contour = vtkSmartPointer<vtkContourFilter>::New();
//    p_contour->SetInput(p_image);
//    p_contour->SetValue(0, 255);
//    p_contour->Update();
//
//    // Connect points in moore neighourhood, do not revist points
//    vtkSmartPointer<vtkPolyData> p_polydata = vtkSmartPointer<vtkPolyData>::New();
//    p_polydata->SetPoints(contour->GetOutput()->GetPoints());
//
//    vtkSmartPointer<vtkKdTreePointLocator> p_locator = vtkSmartPointer<vtkKdTreePointLocator>::New();
//    p_locator->SetDataSet(p_polydata);
//    p_locator->BuildLocator();
//
//    double spacing_x = p_image->GetOutput()->GetSpacing()[0];
//    double spacing_y = p_image->GetOutput()->GetSpacing()[1];
//    double diag_spacing = sqrt(spacing_x*spacing_x + spacing_y*spacing_y);
//
//    vtkSmartPointer<vtkCellArray> p_lines = vtkSmartPointer<vtkCellArray>::New();
//    for(unsigned idx=0; idx<p_polydata->GetNumberOfPoints(); idx++)
//    {
//        vtkSmartPointer<vtkIdList> p_id_list = vtkSmartPointer<vtkIdList>::New();
//        p_locator->FindPointsWithinRadius(diag_spacing + 1.e-6, polydata->GetPoint(idx), p_id_list);
//        for(unsigned jdx=0; jdx<p_id_list->GetNumberOfIds(); jdx++)
//        {
//            if(p_id_list->GetId(jdx) > idx)
//            {
//                vtkSmartPointer<vtkLine> p_line = vtkSmartPointer<vtkLine>::New();
//                p_line->GetPointIds()->InsertId(0, idx);
//                p_line->GetPointIds()->InsertId(1, p_id_list->GetId(jdx));
//                p_lines.InsertNextCell(line);
//            }
//        }
//    }
//    p_polydata->SetLines(lines);
//
//    // So far, this leaves small triangles/kites. Remove them by marking any edges that are 1) diagonal and 2) both points have connectivity 2
//    std::vector<c_vector<double, 2> > points;
//    std::vector<std::vector<unsigned> > edges;
//
//    for(unsigned idx=0; idx<p_polydata->GetNumberOfPoints(); idx++)
//    {
//        c_vector<double, 2> loc;
//        loc[0] = p_polydata->GetPoints(idx)[0];
//        loc[1] = p_polydata->GetPoints(idx)[1];
//        points.push_back(loc);
//    }
//    p_polydata->GetLines()->InitTraversal();
//    vtkSmartPointer<vtkIdList> p_id_list = vtkSmartPointer<vtkIdList>::New();
//    for(unsigned idx=0; idx< p_polydata->GetNumberOfLines(); idx++)
//    {
//        p_polydata->GetLines()->GetNextCell(p_id_list);
//        std::vector<unsigned> point_indices;
//        for(unsigned jdx=0; jdx<p_id_list.GetNumberOfIds(); jdx++)
//        {
//            unsigned seg_id = p_id_list.GetId(jdx);
//            point_indices.push_back(seg_id);
//        }
//        edges.push_back(point_indices);
//    }
//
//    std::vector<std::vector<unsigned> > connectivity(p_polydata->GetNumberOfPoints());
//    for(unsigned idx=0; idx<edges.size(); idx++)
//    {
//        connectivity[edges[idx][0]].push_back(idx);
//        connectivity[edges[idx][1]].push_back(idx);
//    }
//
//    vtkSmartPointer<vtkDoubleArray> p_double_array = vtkSmartPointer<vtkDoubleArray>::New();
//    p_double_array->SetName("ToBeRemoved");
//    p_double_array->.SetNumberOfTuples(edges.size());
//
//    for(unsigned idx=0; idx<edges.size(); idx++)
//    {
//        p_double_array->SetTuple1(idx, 0.0);
//        c_vector<double, 2> loc1 = points[edges[idx][0]];
//        c_vector<double, 2> loc2 = points[edges[idx][1]];
//        double length = sqrt((loc1[0]-loc2[0])* (loc1[0]-loc2[0]) + (loc1[1]-loc2[1])* (loc1[1]-loc2[1]));
//
//        if(length >=diag_spacing - 1.e-6)
//        {
//            // Special case to make sure that pruning kites doesn't disconnect regions
//            if(connectivity[edges[idx][0]].size()>2 and connectivity[edges[idx][1]].size()>2)
//            {
//                unsigned ptid1 = edges[idx][0];
//                unsigned ptid2 = edges[idx][1];
//                // collect all the end point ids. Two of them, not including the end points of this edge
//                // need to be the same to make sure the graph doesn't break due to pruning.
//
//                std::vector<unsigned> end_point_ids;
//                for(unsigned jdx=0; jdx<connectivity[ptid1].size(); jdx++)
//                {
//                    for(unsigned kdx=0; kdx<edges[connectivity[jdx]].size(); kdx++)
//                    {
////                        end_point_ids
//                    }
//                }
//            }
//        }
//    }


//    for idx, eachEdge in enumerate(edges):
//        data.SetTuple1(idx, 0.0)
//        length = np.linalg.norm(np.array(points[eachEdge[0]])- np.array(points[eachEdge[1]]))
//        if length >= diag_spacing - 1.e-6:
//            if len(connectivity[eachEdge[0]]) > 2 and len(connectivity[eachEdge[1]])>2:
//                p1 = eachEdge[0]
//                p2 = eachEdge[1]
//
//                # collect all the end point ids. Two of them, not including the end points of this edge
//                # need to be the same to make sure the graph doesn't break due to pruning.
//                end_point_ids = []
//                for eachEdge in connectivity[p1]:
//                    for eachPoint in edges[eachEdge]:
//                        end_point_ids.append(eachPoint)
//                for eachEdge in connectivity[p2]:
//                    for eachPoint in edges[eachEdge]:
//                        end_point_ids.append(eachPoint)
//                end_point_ids = [value for value in end_point_ids if value != p1]
//                end_point_ids = [value for value in end_point_ids if value != p2]
//                counter=collections.Counter(end_point_ids)
//                if(counter.most_common(1)[0][1])>1:
//                    data.SetTuple1(idx, 1.0)

//    line_boundary.GetCellData().SetScalars(data)

//    threshold = vtk.vtkThreshold()
//    threshold.SetInput(line_boundary)
//    threshold.SetInputArrayToProcess(0, 0, 0, "vtkDataObject::FIELD_ASSOCIATION_CELLS", "ToBeRemoved")
//    threshold.ThresholdByLower(0.1)
//    threshold.Update()
//
//    polysurface = vtk.vtkGeometryFilter()
//    polysurface.SetInputConnection(threshold.GetOutputPort())
//    polysurface.Update()
//
//    # Join individual lines into polylines
//    stripper = vtk.vtkStripper()
//    stripper.SetInput(polysurface.GetOutput())
//    stripper.Update()
//    stripped = stripper.GetOutput()
//    #
//    # # Label each line
//    data = vtk.vtkDoubleArray()
//    data.SetNumberOfTuples(stripped.GetNumberOfLines())
//    for idx in range(stripped.GetNumberOfLines()):
//        data.SetTuple1(idx, idx)
//    stripped.GetCellData().SetScalars(data)
//    chaste.utility.readwrite.write(stripped, "/home/grogan/with_lines_stripped.vtp")


}
#endif /*CHASTE_ANGIOGENESIS_EXTENDED*/
