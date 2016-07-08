import vtk
import numpy as np
from meshpy.triangle import MeshInfo, build
#import dolfin as df
import numpy as np
import chaste.utility.bases as bases
import chaste.mesh
import chaste.mesh.converters
import chaste.population.vessel
import chaste.image

# class VesselNetworkMesher1d(bases.SimpleIOBase):
#     
#     def __init__(self):
#         super(VesselNetworkMesher1d, self).__init__()
#         
#         self.mesh_size = 10.0
#         self.su_subdivisions = 10
#         
#     def update(self):
#         
#         # Get the vtk representation
#         writer = chaste.population.vessel.VtkVesselNetworkWriter()
#         writer.SetVesselNetwork(self.input)
#         vtk_rep = writer.GetOutput()
#         
#         # subdivide
#         spline = vtk.vtkSplineFilter();
#         spline.SetInput(vtk_rep)
# #        spline.SetLength(self.mesh_size*10.0)
#         spline.SetNumberOfSubdivisions(10)
#         
#         spline.Update();
#         
#         triangle = vtk.vtkTriangleFilter()
#         triangle.SetInput(spline.GetOutput())
#         triangle.Update()
#         
#         clean = vtk.vtkCleanPolyData()
#         clean.SetInput(triangle.GetOutput())
#         clean.Update()
#         
#         vtk_rep = clean.GetOutput();
#         
#         # convert vtk to MeshPy Tri format
#         converter = converters.VtkToTriMesh()
#         converter.input = vtk_rep
#         converter.update()
#         points, edges = converter.output
#         
#         editor = df.MeshEditor()
#         centre_mesh = df.Mesh()
#         
#         editor.open(centre_mesh, 1, 2)
#         editor.init_vertices(len(points))
#         editor.init_cells(len(edges))
# 
#         for idx, eachPoint in enumerate(points):
#             editor.add_vertex(idx, np.array((eachPoint[0], eachPoint[1]), dtype=np.float))
#             
#         for idx, eachCell in enumerate(edges):
#             editor.add_cell(idx, np.array((eachCell[0], eachCell[1]), dtype=np.uintp))
# 
#         editor.close()
#         centre_mesh.init()
#         centre_mesh.order()
#         
#         file = df.File("/home/grogan/ABME16_Work/test/TestSingleVesselMeshing/mesh_1d.xml")
#         file << centre_mesh
#         return centre_mesh
    
class NetworkToMesh(bases.SimpleIOBase):
     
    def __init__(self):
        super(NetworkToMesh, self).__init__()
         
        self.mesh_size = 10.0
        self.su_subdivisions = 10
         
    def update(self):
         
        image_converter = chaste.image.NetworkToImage()
        
        
        image_converter.SetNetwork(self.input)
        image_converter.SetGridSpacing(2.0)
        image_converter.SetPaddingFactors(0.1, 0.1, 0.0)
        image_converter.SetImageDimension(2)
        image_converter.Update()
        image = image_converter.GetOutput()
           
        marching_squares = vtk.vtkMarchingSquares()
        marching_squares.SetInput(image)
        marching_squares.SetValue(0, 1)
        marching_squares.Update()
        
        cleaner = vtk.vtkCleanPolyData()
        cleaner.SetInput(marching_squares.GetOutput())
        cleaner.Update()
        
        stripper = vtk.vtkStripper()
        stripper.SetInput(cleaner.GetOutput())
        stripper.Update()
        
        # Downsample and smooth the polyline
        spline = vtk.vtkSplineFilter()
        spline.SetLength(10.0)
        spline.SetSubdivideToLength()
        spline.SetInput(stripper.GetOutput())
        spline.Update()
        
        triangle = vtk.vtkTriangleFilter()
        triangle.SetInput(spline.GetOutput())
        triangle.Update()
        
        cleaned = triangle.GetOutput()
        # Want flat ends on input and output network nodes. Do this after smoothing
        for eachNode in self.input.GetNodes():
            if eachNode.GetFlowProperties().IsInputNode() or eachNode.GetFlowProperties().IsOutputNode():
                
                # Set up a generic clipping box
                radius = eachNode.GetRadius()
                box = vtk.vtkBox()
                box.SetBounds(-1.1*radius, 0.0, -1.1*radius, 1.1*radius, - 1.1*radius, 1.1*radius)
                
                # Align the box with the vessel
                loc = np.array(eachNode.GetLocation())
                tangent = np.array(eachNode.GetSegment(0).GetOppositeNode(eachNode).GetLocation()) - loc
                tangent /= np.linalg.norm(tangent)
                
                box_axis = np.array((1.0, 0.0, 0.0))
                rotation_angle = np.rad2deg(np.arccos(np.dot(box_axis, tangent)))  
                rotation_axis = np.cross(box_axis, tangent)
                

                transform = vtk.vtkTransform()
                if(np.abs(np.dot(box_axis, tangent))<1.0-1.e-6):
                    transform.RotateWXYZ(-rotation_angle, rotation_axis[0], rotation_axis[1], rotation_axis[2])
                else:
                    transform.RotateWXYZ(-rotation_angle, 0.0, 0.0, 1.0)
                translation_vector = -loc
                print translation_vector
                transform.Translate(translation_vector[0],translation_vector[1], translation_vector[2])
                box.SetTransform(transform)
                
                clipper = vtk.vtkClipPolyData()
                clipper.SetInput(cleaned)
                clipper.SetClipFunction(box)
                clipper.Update()
                clipped = clipper.GetOutput()
                
                # Assuming we have two points with connectivity 1, join them
                edge_ids = []
                for idx in range(clipped.GetNumberOfPoints()):
                    cell_list = vtk.vtkIdList()
                    clipped.GetPointCells(idx, cell_list)
                    if(cell_list.GetNumberOfIds()==1):
                        edge_ids.append(idx)
                         
                line = vtk.vtkLine()
                line.GetPointIds().SetId(0, edge_ids[0])
                line.GetPointIds().SetId(1, edge_ids[1])
                clipped.GetLines().InsertNextCell(line)
                clipped.Update()
                
                cleaned = vtk.vtkPolyData()
                cleaned.DeepCopy(clipped)
                
#        return cleaned
                
        tri_input_converter = chaste.mesh.converters.VtkToTriMesh()
        tri_input_converter.input = cleaned
        points, edges = tri_input_converter.update()
         
        # First generate a coarse mesh and use it to probe for holes
        mesh_info = MeshInfo()
        mesh_info.set_points(points)
        mesh_info.set_facets(edges)
         
        # Allow for two different regions and holes
        data = build(mesh_info)
        points= data.points
        edges = data.elements
         
        vtk_points = vtk.vtkPoints()
        for eachEdge in edges:
            centx = (points[eachEdge[0]][0] + points[eachEdge[1]][0] + points[eachEdge[2]][0])/3.0
            centy = (points[eachEdge[0]][1] + points[eachEdge[1]][1] + points[eachEdge[2]][1])/3.0
            vtk_points.InsertNextPoint(centx, centy, 0.0)
             
        # Get the values of the image data at the points
        temp_polydata = vtk.vtkPolyData()
        temp_polydata.SetPoints(vtk_points)
        probe = vtk.vtkProbeFilter()
        probe.SetInput(temp_polydata)
        probe.SetSource(image)
        probe.Update()
        results = probe.GetOutput().GetPointData().GetScalars()
         
        hole_locations = []
        for idx in range(results.GetNumberOfTuples()):
            if(results.GetTuple1(idx)==0):
                loc = vtk_points.GetPoint(idx)
                hole_locations.append([loc[0], loc[1]])
         
        if(len(hole_locations)>0):
            mesh_info.holes.resize(len(hole_locations))
            for idx, eachHole in enumerate(hole_locations):
                mesh_info.holes[idx] = [eachHole[0], eachHole[1]]
                 
        # Remesh with holes
        data = build(mesh_info, max_volume = self.mesh_size)
        tri_to_vtk = chaste.mesh.converters.TriMeshToVtkUnstructured()
        tri_to_vtk.input = [data.points, data.elements]
        mesh = tri_to_vtk.update()
        return mesh
    
class NetworkToMesh3d(bases.SimpleIOBase):
     
    def __init__(self):
        super(NetworkToMesh3d, self).__init__()
         
        self.mesh_size = 10.0
        self.su_subdivisions = 10
         
    def update(self):
         
        # Get the image representation, get the edges, smooth, mesh 2d
        image_converter = chaste.image.NetworkToImage()
        image_converter.SetNetwork(self.input)
        image_converter.SetGridSpacing(2.0)
        image_converter.SetPaddingFactors(0.1, 0.1, 0.1)
        image_converter.SetImageDimension(3)
        image_converter.Update()
        image = image_converter.GetOutput()
           
        marching_cubes = vtk.vtkMarchingCubes()
        marching_cubes.SetInput(image)
        marching_cubes.SetValue(0, 1)
        marching_cubes.Update()
        
        cleaner = vtk.vtkCleanPolyData()
        cleaner.SetInput(marching_cubes.GetOutput())
        cleaner.Update()
        
        smoother = vtk.vtkWindowedSincPolyDataFilter()
        smoother.SetInput(cleaner.GetOutput())
        smoother.SetNumberOfIterations(100.0)
        smoother.BoundarySmoothingOn()
        smoother.FeatureEdgeSmoothingOn()        
        smoother.SetFeatureAngle(130)
        smoother.SetPassBand(0.05) 
        smoother.NonManifoldSmoothingOn() 
        smoother.NormalizeCoordinatesOn() 
        smoother.Update()
        
        cleaned = smoother.GetOutput()
        # Want flat ends on input and output network nodes. Do this after smoothing
        for eachNode in self.input.GetNodes():
            if eachNode.GetFlowProperties().IsInputNode() or eachNode.GetFlowProperties().IsOutputNode():
                
                # Set up a generic clipping box
                radius = eachNode.GetRadius()
                box = vtk.vtkBox()
                box.SetBounds(-1.1*radius, 0.0, -1.1*radius, 1.1*radius, - 1.1*radius, 1.1*radius)
                
                # Align the box with the vessel
                loc = np.array(eachNode.GetLocation())
                tangent = np.array(eachNode.GetSegment(0).GetOppositeNode(eachNode).GetLocation()) - loc
                tangent /= np.linalg.norm(tangent)
                
                box_axis = np.array((1.0, 0.0, 0.0))
                rotation_angle = np.rad2deg(np.arccos(np.dot(box_axis, tangent)))  
                rotation_axis = np.cross(box_axis, tangent)
                
                transform = vtk.vtkTransform()
                if(np.abs(np.dot(box_axis, tangent))<1.0-1.e-6):
                    transform.RotateWXYZ(-rotation_angle, rotation_axis[0], rotation_axis[1], rotation_axis[2])
                else:
                    transform.RotateWXYZ(-rotation_angle, 0.0, 0.0, 1.0)
                translation_vector = -loc
                print translation_vector
                transform.Translate(translation_vector[0],translation_vector[1], translation_vector[2])
                box.SetTransform(transform)
                
                clipper = vtk.vtkClipPolyData()
                clipper.SetInput(cleaned)
                clipper.SetClipFunction(box)
                clipper.Update()
                
#                 bound_edges = vtk.vtkFeatureEdges()
#                 bound_edges.SetInput(clipper.GetOutput())
#                 bound_edges.BoundaryEdgesOn()
#                 bound_edges.FeatureEdgesOff()
#                 bound_edges.NonManifoldEdgesOff()
#                 bound_edges.ManifoldEdgesOff()
#                 
#                 strips = vtk.vtkStripper()
#                 strips.SetInput(bound_edges.GetOutput())
#                 strips.Update()
#                 
#                 bound_polys = vtk.vtkPolyData()
#                 bound_polys.SetPoints(strips.GetOutput().GetPoints())
#                 bound_polys.SetPolys(strips.GetOutput().GetLines())
#                 
#                 append = vtk.vtkAppendPolyData()
#                 append.AddInput(clipper.GetOutput())
#                 append.AddInput(bound_polys)
#                 append.Update()
#                 
#                 triangle = vtk.vtkTriangleFilter()
#                 triangle.SetInput(append.GetOutput())
#                 triangle.Update()
#                 
#                 clean= vtk.vtkCleanPolyData()
#                 clean.SetInput(triangle.GetOutput())
#                 clean.Update()
                
                cleaned.DeepCopy(clipper.GetOutput())        # Clip ends
        
              
        return cleaned
        