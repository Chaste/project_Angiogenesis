"""
Single vessel network
"""

import numpy as np
import vtk
from vtk.util import numpy_support
import chaste.population.vessel

class PartToGeometry():
    
    def __init__(self, part):
        
        self.Base = __import__("FreeCAD.Base")
        self.Part = __import__("Part")

        self.part = part
        
    def generate(self):
        
        polygons = self.part.GetPolygons()
        faces = []
        for eachPolygon in polygons:
            locs = []
            verts = eachPolygon.get_vertices()
            for eachVert in verts:
                locs.append(tuple(eachVert.get_location()))
            locs.append(verts[0].get_location())
            poly = self.Part.makePolygon(locs)
            
            faces.append(self.Part.Face(poly))
            
        shell = self.Part.makeShell(faces)
        
        return shell

class NetworkToGeometry():
    
    def __init__(self):
        
        self.Base = __import__("FreeCAD.Base")
        self.Part = __import__("Part")
        
        self.network = None
        self.dimension = 3
        self.geometry = None
        
    def set_network(self, network):
        self.network = network
        
    def set_dimension(self, dimension):
        self.dimension = dimension
        
    def update(self):
        
        # make a cylinder for each vessel
        vessels = self.network.vessels
        
        for eachVessel in vessels:
            
            start = eachVessel.start_node.GetLocationVector()
            end = eachVessel.end_node.GetLocationVector()
        
            cylinder = self.Part.makeCylinder(eachVessel.start_node.radius, eachVessel.length, 
                                         self.Base.Vector(tuple(start)), self.Base.Vector(tuple(end - start)))
            
            cylinder = self.Part.makeShell(cylinder.Faces)
            
        if self.dimension == 2:
            
            clipping_plane = self.Part.makePlane(1.e6, 1.e6, self.Base.Vector(-5.e5,  -5.e5, 0.0), self.Base.Vector(0,0,1))
            section = cylinder.section(clipping_plane)
            
            polyData = vtk.vtkPolyData() 
            points = vtk.vtkPoints()
            lines = vtk.vtkCellArray()
            
            for eachEdge in section.Edges:
                line = vtk.vtkLine()
                pointId1 = points.InsertNextPoint(eachEdge.Vertexes[0].Point)
                pointId2 = points.InsertNextPoint(eachEdge.Vertexes[1].Point)
                line.GetPointIds().InsertId(0, pointId1)
                line.GetPointIds().InsertId(1, pointId2)
                lines.InsertNextCell(line)
                
            polyData.SetPoints(points)
            polyData.SetLines(lines)   
            
            clean = vtk.vtkCleanPolyData()
            clean.SetInput(polyData)
            clean.Update()
        
            cylinder = clean.GetOutput()
            
        self.geometry = cylinder
            
        return self.geometry
    
    def get_output(self):
        return self.geometry
    
class NetworkToVtk():
    
    def __init__(self):
        
        self.network = None
        
    def set_network(self, network):
        self.network = network
        
    def generate(self):
        polyData = vtk.vtkPolyData() 
        points = vtk.vtkPoints()
        lines = vtk.vtkCellArray()
        info = vtk.vtkFloatArray()
        
        info.SetNumberOfComponents(1)
        info.SetName("Radius")
        
        for eachVessel in self.network.vessels:
            line = vtk.vtkLine()
            avRadius = 0
            for idx, eachSegment in enumerate(eachVessel.vesselSegments):
                pointId = points.InsertNextPoint(eachSegment.startNode.location)
                line.GetPointIds().InsertId(idx, pointId)
                avRadius = avRadius + eachSegment.radius
                if idx == len(eachVessel.vesselSegments) - 1:
                    pointId = points.InsertNextPoint(eachSegment.endNode.location)
                    line.GetPointIds().InsertId(idx + 1, pointId)
            avRadius = avRadius / len(eachVessel.vesselSegments)
            lines.InsertNextCell(line)
            info.InsertNextTupleValue([avRadius])   
                
        polyData.SetPoints(points)
        polyData.SetLines(lines)
        polyData.GetCellData().SetScalars(info)
        
        return polyData
    
class TriToVtk():
    
    def __init__(self, mesh):
        
        self.mesh = mesh
        
    def generate(self):
        
        points = self.mesh[0]
        edges = self.mesh[1]
        polyData = vtk.vtkPolyData() 
        new_vtk_points = vtk.vtkPoints()
        new_vtk_lines = vtk.vtkCellArray()
        
        for eachEdge in edges:
            line = vtk.vtkLine()
            loc1 = (points[eachEdge[0]][0], points[eachEdge[0]][1], 0.0)
            loc2 = (points[eachEdge[1]][0], points[eachEdge[1]][1], 0.0)
            pointId1 = new_vtk_points.InsertNextPoint(loc1)
            pointId2 = new_vtk_points.InsertNextPoint(loc2)
            line.GetPointIds().InsertId(0, pointId1)
            line.GetPointIds().InsertId(1, pointId2)
            new_vtk_lines.InsertNextCell(line)
                
        polyData.SetPoints(new_vtk_points)
        polyData.SetLines(new_vtk_lines)
        
        return polyData
    
class VtkToTri():
    
    def __init__(self):
        
        self.polydata = None
        self.points = None
        self.edges = None
        
    def generate(self, polydata):
        
        self.polydata = polydata
        vtk_numPoints = self.polydata.GetNumberOfPoints()    
        vtk_points = self.polydata.GetPoints()  
    
        self.points= []
        self.edges = []
        for i in range(vtk_numPoints):
            self.points.append([vtk_points.GetPoint(i)[0], vtk_points.GetPoint(i)[1]])
               
        numCells = self.polydata.GetNumberOfLines()  
        cellArray = self.polydata.GetLines()
        cellArray.InitTraversal()
        segList = vtk.vtkIdList()
        
        for i in range(numCells): 
            cellArray.GetNextCell(segList)
            point_indices = []
            for j in range(0, segList.GetNumberOfIds()):
                seg_id = segList.GetId(j)
                point_indices.append(int(seg_id))
            self.edges.append(point_indices)
            
        return self.points, self.edges
        
class VtkToNetwork():
    
    def __init__(self, polydata):
        
        self.polydata = polydata
        
    def generate(self):
        
        vtk_numPoints = self.polydata.GetNumberOfPoints()    
        vtk_points = self.polydata.GetPoints()  
        
        nodes = [] 
        for i in range(vtk_numPoints):
            nodes.append(chaste.population.vessel.VascularNode(vtk_points.GetPoint(i)))
        
        numCells = self.polydata.GetNumberOfLines()  
        cellArray = self.polydata.GetLines()
        cellArray.InitTraversal()
        segList = vtk.vtkIdList()
        
        network = chaste.population.vessel.VascularNetwork()
        for i in range(numCells): 
            cellArray.GetNextCell(segList)
            point_indices = []
            for j in range(0, segList.GetNumberOfIds()):
                seg_id = segList.GetId(j)
                point_indices.append(int(seg_id))
            vessel = chaste.population.vessel.Vessel([nodes[point_indices[0]], nodes[point_indices[1]]])
            network.addVessel(vessel)   
        network.UpdateAll(False)
        
        return network
    
class CellConverter:
    
    def __init__(self):
        pass   
    
    def convert(self, cell_population):
        
        polyData = vtk.vtkPolyData() 
        points = vtk.vtkPoints()

        for eachCell in cell_population.get_cells():
            points.InsertNextPoint(eachCell.get_location())  
        polyData.SetPoints(points)
        
        return polyData
    
class PointConverter:
    
    def __init__(self):
        pass   
    
    def convert(self, points):
        
        polyData = vtk.vtkPolyData() 
        vtk_points = vtk.vtkPoints()
        for eachPoint in points:
            vtk_points.InsertNextPoint(eachPoint.get_location())    
        polyData.SetPoints(vtk_points)
        
        return polyData
    
class MeshConverter():
    
    def __init__(self):
        pass
    
    def convert(self, mesh, dimension = 2):
        
        grid = vtk.vtkUnstructuredGrid()
        
        # Add VTK points corresponding to mesh nodes
        points = vtk.vtkPoints()
        locations = mesh.get_node_locations()
        points.SetNumberOfPoints(len(locations))
        for idx, eachLocation in enumerate(locations):
            args = [idx] + list(eachLocation)
            points.InsertPoint(*args)
        grid.SetPoints(points)  

        # Add VTK Tets or Triangles corresponding to mesh elements
        connectivity = mesh.get_connectivity()
        num_elements = len(connectivity)
        grid.Allocate(num_elements, num_elements)
        for idx in range(num_elements): 
            if dimension == 3:  
                vtkElement = vtk.vtkTetra()
            else:
                vtkElement = vtk.vtkTriangle()

            num_nodes = len(connectivity[idx])
            for jdx in range(num_nodes):  
                vtkElement.GetPointIds().SetId(jdx, connectivity[idx][jdx])                       
            grid.InsertNextCell(vtkElement.GetCellType(), vtkElement.GetPointIds())
        return grid
    
def vtkpolygon_to_lines(surface):
    
    feature_edges = vtk.vtkFeatureEdges()
    feature_edges.SetInput(surface)  
    feature_edges.Update()   
        
    clean = vtk.vtkCleanPolyData()
    clean.SetInputConnection(feature_edges.GetOutputPort())
    clean.Update()
    
    triangle2 = vtk.vtkTriangleFilter()
    triangle2.SetInputConnection(clean.GetOutputPort())
    triangle2.Update()

    return triangle2.GetOutput()

def vtk_network_lines_to_polylines(boundary):
        
    numPoints = boundary.GetNumberOfPoints()    
    vtkpoints = boundary.GetPoints()  
    connectivity = []
    
    points = [] 
    for i in range(numPoints):
        points.append(vtkpoints.GetPoint(i))
        connectivity.append([])
        
    numCells = boundary.GetNumberOfLines()  
    cellArray = boundary.GetLines()
    cellArray.InitTraversal()
    segList = vtk.vtkIdList()
        
    edges = []
    for i in range(numCells): 
        cellArray.GetNextCell(segList)
        point_indices = []
        for j in range(0, segList.GetNumberOfIds()):
            seg_id = segList.GetId(j)
            point_indices.append(int(seg_id))
        edges.append((point_indices[0], point_indices[1]))
        connectivity[point_indices[0]].append(i)
        connectivity[point_indices[1]].append(i)
        
    branches = []
    point_visited = np.zeros(len(points))
    
    for idx, eachPoint in enumerate(points):
        num_segs = len(connectivity[idx])

        if num_segs ==3 or num_segs==1:
            point_visited[idx] = 1
            for eachEdgeId in connectivity[idx]:
                branch_points = [idx]
                current_point_id = idx
                previous_point_id = idx
                found_branch = False
                while not found_branch:
                      
                    # get the next node   
                    current_num_segs = len(connectivity[current_point_id])
                    if current_num_segs == 2:
                        edge_1 = edges[connectivity[current_point_id][0]]
                        edge_2 = edges[connectivity[current_point_id][1]]
                        if edge_1[0] == current_point_id:
                            opp1 = edge_1[1]
                        else:
                            opp1 = edge_1[0]
                        if opp1 != previous_point_id:
                            next_edge = edge_1
                        else:
                            next_edge = edge_2
                    else:
                        next_edge = edges[eachEdgeId]
                              
                    if next_edge[0]== current_point_id:
                        next_point_id = next_edge[1]  
                    else:
                        next_point_id = next_edge[0]     
                              
                    if point_visited[next_point_id] == 1 and len(connectivity[next_point_id]) == 2:
                        break
                      
                    point_visited[next_point_id] = 1
                    branch_points.append(next_point_id)
                    
                    if len(connectivity[next_point_id]) == 2:
                        previous_point_id = current_point_id
                        current_point_id = next_point_id
                    else:
                        found_branch = True
                        break
                          
                if len(branch_points) > 1:
                    branches.append(branch_points)
                    
    new_points = vtk.vtkPoints()
    new_points.SetNumberOfPoints(len(points))
    for idx, eachPoint in enumerate(points):
        new_points.SetPoint(idx, points[idx][0], points[idx][1], 0.0)
        
    lines = vtk.vtkCellArray()
    for eachBranch in branches:
        lines.InsertNextCell(len(eachBranch))
        for eachPoint in eachBranch:
            lines.InsertCellPoint(eachPoint)
            
    polygon = vtk.vtkPolyData()
    polygon.SetPoints(new_points)
    polygon.SetLines(lines)
    
    return polygon  
    
class ImageDataToNumpy():
    
    def __init__(self):
        pass
    
    def convert(self, data):
        spacing = data.GetSpacing()
        dims = data.GetDimensions()
        origin = data.GetOrigin()
        vtk_data = data.GetPointData().GetScalars()
         
        numpy_data = numpy_support.vtk_to_numpy(vtk_data)
        numpy_data = numpy_data.reshape(dims[0], dims[1], dims[2])
        numpy_data = numpy_data.transpose(2,1,0)
        
        return origin, dims, spacing, numpy_data
    