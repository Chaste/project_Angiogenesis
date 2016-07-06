import numpy as np
import collections
import vtk
import chaste.utility.readwrite
import chaste.mesh.converters

def network_from_image(image_file, work_dir, prune_level = 2):

    image = chaste.utility.readwrite.read(image_file)
    
    
    # get the edges
    threshold = vtk.vtkThreshold()
    threshold.SetInput(image)
    threshold.ThresholdByUpper(0.1)
    threshold.Update()
    
    polysurface = vtk.vtkGeometryFilter()
    polysurface.SetInputConnection(threshold.GetOutputPort())
    polysurface.Update()
    
    features = vtk.vtkFeatureEdges()
    features.SetInputConnection(polysurface.GetOutputPort())
    chaste.utility.readwrite.write(features.GetOutput(), work_dir + "edges.vtp")
    
    # Get the skeleton
    skeleton =vtk.vtkImageSkeleton2D()
    skeleton.SetInput(image)
    skeleton.SetNumberOfIterations(300)
    skeleton.SetPrune(prune_level)
    skeleton.ReleaseDataFlagOff()
    skeleton.Update()
    chaste.utility.readwrite.write(skeleton.GetOutput(), work_dir + "skeleton.vti")
    
    # Use a contour filter to get points on the skeleton
    contour = vtk.vtkContourFilter()
    contour.SetInput(skeleton.GetOutput())
    contour.SetValue(0, 255)
    contour.Update()
    
    # Connect points in moore neighourhood, do not revist points
    polydata = vtk.vtkPolyData()
    polydata.SetPoints(contour.GetOutput().GetPoints())
    locator = vtk.vtkKdTreePointLocator()
    locator.SetDataSet(polydata)
    locator.BuildLocator()
    
    spacing_xy = np.array((skeleton.GetOutput().GetSpacing()[0], skeleton.GetOutput().GetSpacing()[1]))
    diag_spacing  = np.linalg.norm(spacing_xy)
    lines = vtk.vtkCellArray()
    for idx in range(polydata.GetNumberOfPoints()):
        idList = vtk.vtkIdList()
        locator.FindPointsWithinRadius(diag_spacing + 1.e-6, polydata.GetPoint(idx), idList)
        for jdx in range(idList.GetNumberOfIds()):
            if idList.GetId(jdx) > idx:
                line = vtk.vtkLine()
                line.GetPointIds().InsertId(0, idx)
                line.GetPointIds().InsertId(1, idList.GetId(jdx))
                lines.InsertNextCell(line)
    polydata.SetLines(lines)

    # So far, this leaves small triangles. Remove them by marking any edges that are 1) diagonal and 2) both points have connectivity 2
    line_boundary = polydata
    converter = chaste.mesh.converters.VtkToTriMesh()
    converter.input = line_boundary
    points, edges = converter.update()
    
    connectivity = []
    for idx in range(len(points)):
        connectivity.append([])
    for idx, eachEdge in enumerate(edges):
        connectivity[eachEdge[0]].append(idx)
        connectivity[eachEdge[1]].append(idx)
          
    data = vtk.vtkDoubleArray()
    data.SetName("ToBeRemoved")
    data.SetNumberOfTuples(len(edges))
    
    for idx, eachEdge in enumerate(edges):
        data.SetTuple1(idx, 0.0)
        length = np.linalg.norm(np.array(points[eachEdge[0]])- np.array(points[eachEdge[1]]))
        if length >= diag_spacing - 1.e-6:
            if len(connectivity[eachEdge[0]]) > 2 and len(connectivity[eachEdge[1]])>2:
                p1 = eachEdge[0]
                p2 = eachEdge[1]
                
                # collect all the end point ids. Two of them, not including the end points of this edge
                # need to be the same to make sure the graph doesn't break due to pruning.
                end_point_ids = []
                for eachEdge in connectivity[p1]:
                    for eachPoint in edges[eachEdge]:
                        end_point_ids.append(eachPoint)
                for eachEdge in connectivity[p2]:
                    for eachPoint in edges[eachEdge]:
                        end_point_ids.append(eachPoint)
                end_point_ids = [value for value in end_point_ids if value != p1]
                end_point_ids = [value for value in end_point_ids if value != p2]
                counter=collections.Counter(end_point_ids)
                if(counter.most_common(1)[0][1])>1:          
                    data.SetTuple1(idx, 1.0)
            
    line_boundary.GetCellData().SetScalars(data)
    
    threshold = vtk.vtkThreshold()
    threshold.SetInput(line_boundary)
    threshold.SetInputArrayToProcess(0, 0, 0, "vtkDataObject::FIELD_ASSOCIATION_CELLS", "ToBeRemoved")
    threshold.ThresholdByLower(0.1)
    threshold.Update()
            
    polysurface = vtk.vtkGeometryFilter()
    polysurface.SetInputConnection(threshold.GetOutputPort())
    polysurface.Update()
    
    # Join individual lines into polylines
    stripper = vtk.vtkStripper()
    stripper.SetInput(polysurface.GetOutput())
    stripper.Update()
    stripped = stripper.GetOutput()
    # 
    # # Label each line
    data = vtk.vtkDoubleArray()
    data.SetNumberOfTuples(stripped.GetNumberOfLines())
    for idx in range(stripped.GetNumberOfLines()):
        data.SetTuple1(idx, idx)
    stripped.GetCellData().SetScalars(data)
    
    extrude = vtk.vtkLinearExtrusionFilter()
    extrude.SetInputConnection( features.GetOutputPort());
    extrude.SetExtrusionTypeToNormalExtrusion();
    extrude.SetVector(0, 0, 1 );
    extrude.SetScaleFactor (1.0);
    extrude.Update();
     
    triangleFilter = vtk.vtkTriangleFilter();
    triangleFilter.SetInputConnection(extrude.GetOutputPort());
    triangleFilter.Update();
    
    extrude2 = vtk.vtkLinearExtrusionFilter()
    extrude2.SetInput(stripped);
    extrude2.SetExtrusionTypeToNormalExtrusion();
    extrude2.SetVector(0, 0, 1 );
    extrude2.SetScaleFactor (1.0);
    extrude2.Update();
     
    triangleFilter2 = vtk.vtkTriangleFilter();
    triangleFilter2.SetInputConnection(extrude2.GetOutputPort());
    triangleFilter2.Update();
    
    distance = vtk.vtkDistancePolyDataFilter()
    distance.SetInput(0, triangleFilter2.GetOutput());
    distance.SetInput(1, triangleFilter.GetOutput());
    distance.Update()
    
    plane = vtk.vtkPlane()
    plane.SetOrigin(0.0, 0, 0.5)
    plane.SetNormal(0, 0, 1)
    cutter = vtk.vtkCutter()
    cutter.SetInputConnection(distance.GetOutputPort())
    cutter.SetCutFunction(plane)
    cutter.Update()
    
    stripper = vtk.vtkStripper()
    stripper.SetInput(cutter.GetOutput())
    stripper.Update()
    stripped = stripper.GetOutput()
    
    p_spline = vtk.vtkSplineFilter();
    p_spline.SetInput(stripped);
    p_spline.SetLength(50.0)
    p_spline.Update()
            
    chaste.utility.readwrite.write(p_spline.GetOutput(), work_dir + "distance.vtp")
