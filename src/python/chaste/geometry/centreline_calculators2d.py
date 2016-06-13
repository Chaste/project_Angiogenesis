import numpy as np
import vtk
from vmtk import vmtkscripts, pypes
from scipy.spatial import Voronoi
from scipy.sparse import csr_matrix
import scipy.sparse.csgraph
import chaste.geometry.converters
import chaste.interfaces.vtk_tools.vtk_tools

class Centrelines2d():
    
    def __init__(self):
        self.surface = None
        self.network = None
        self.start_points= None
        
    def set_surface(self, surface):
        
        self.surface = surface

    def update(self):
        
        # Convert from vtk format to tri plc format
        vtk_to_tri = chaste.geometry.converters.VtkToTri()
        points, edges = vtk_to_tri.generate(self.surface)
        
        # Get the voronoi diagram of the input points
        input_points = []
        for eachEdge in edges:
            input_points.append([eachEdge[0], eachEdge[1]])
        vor = Voronoi(points)
        verts = vor.vertices
        ridges = vor.ridge_vertices
        
        # Get the finite ridges
        finite_ridges = []
        for eachRidge in ridges:
            if eachRidge[0] >= 0 and eachRidge[1] >= 0:
                finite_ridges.append(eachRidge)
        finite_ridges
        
        # Only keep ridges which do not intersect the original polygon
        internal_ridges = []
        for eachRidge in finite_ridges:
            
            crosses = False
            for eachEdge in edges:
                p1 = (points[eachEdge[0]][0], points[eachEdge[0]][1], 0.0)
                p2 = (points[eachEdge[1]][0], points[eachEdge[1]][1], 0.0)
                x1 = (verts[eachRidge[0]][0], verts[eachRidge[0]][1], 0.0)
                x2 = (verts[eachRidge[1]][0], verts[eachRidge[1]][1], 0.0)
    
                para1 = vtk.mutable(0)
                para2 = vtk.mutable(0)
                intersect = vtk.vtkLine.Intersection(p1, p2, x1, x2, para1, para2)
                if intersect != 0:
                    crosses = True
                    break
            if not crosses:
                internal_ridges.append(eachRidge)
                
        # Find the largest connected path
        G_dense = np.zeros((len(verts), len(verts)))
        G_verts = []
        for idx in range(len(verts)):
            G_verts.append([])
        
        for eachRidge in internal_ridges:
            G_verts[eachRidge[0]].append(eachRidge[1])
            G_verts[eachRidge[1]].append(eachRidge[0])
            G_dense[eachRidge[0], eachRidge[1]] = 1.0
            G_dense[eachRidge[1], eachRidge[0]] = 1.0
            
        G_sparse = csr_matrix(G_dense)
        num, labels = scipy.sparse.csgraph.connected_components(G_sparse, False)
        label_count = [np.sum(labels == i) for i in range(num)]
        max_label = label_count.index(max(label_count))
        sparse_edge_labels =  np.where(labels == max_label)
        
        sparse_edges = []
        for eachLabel in sparse_edge_labels[0]:
            for eachNeighbour in G_verts[eachLabel]:
                sparse_edges.append([eachLabel, eachNeighbour])
                
        # Use VTK to remove duplicate points
        tri_to_vtk = chaste.geometry.converters.TriToVtk([verts, sparse_edges])
        polyData = tri_to_vtk.generate()
        
        clean = vtk.vtkCleanPolyData()
        clean.SetInput(polyData)
        clean.Update()
        
        # Remove duplicate edges
        centre = clean.GetOutput()
        vtk_to_tri = chaste.geometry.converters.VtkToTri()
        points, edges = vtk_to_tri.generate(centre)
        unique_edges = []

        for idx, eachEdge in enumerate(edges):
            duplicate = False
            for eachInnerEdge in unique_edges:
                if eachInnerEdge[0] == eachEdge[0] and eachInnerEdge[1] == eachEdge[1]:
                    duplicate = True
                    break
                if eachInnerEdge[1] == eachEdge[0] and eachInnerEdge[0] == eachEdge[1]:
                    duplicate = True
                    break 
            if not duplicate:
                unique_edges.append(eachEdge)
        new_vtk_lines = vtk.vtkCellArray()
        for eachEdge in unique_edges:
            line = vtk.vtkLine()
            line.GetPointIds().InsertId(0, eachEdge[0])
            line.GetPointIds().InsertId(1, eachEdge[1])
            new_vtk_lines.InsertNextCell(line)
        
        tidy_edges = vtk.vtkPolyData() 
        tidy_edges.SetPoints(centre.GetPoints())
        tidy_edges.SetLines(new_vtk_lines)
        
        clean = vtk.vtkCleanPolyData()
        clean.SetInput(tidy_edges)
        clean.Update()
        tidy_edges = clean.GetOutput()
        
        vtk_numPoints = tidy_edges.GetNumberOfPoints()    
        vtk_points = tidy_edges.GetPoints() 
        cellArray = tidy_edges.GetLines()
        locator = vtk.vtkKdTreePointLocator()
        locator.SetDataSet(tidy_edges)
        
        for idx, eachLoc in enumerate(self.start_points):
            probe_loc = np.array((eachLoc[0], eachLoc[1], 0.0))
            closest_id = locator.FindClosestPoint(probe_loc)
            vtk_points.InsertNextPoint(probe_loc)
            
            cellArray.InsertNextCell(2)
            cellArray.InsertCellPoint(closest_id)
            cellArray.InsertCellPoint(vtk_numPoints + idx)
        polygon = vtk.vtkPolyData()
        polygon.SetPoints(vtk_points)
        polygon.SetLines(cellArray)
        
        clean = vtk.vtkCleanPolyData()
        clean.SetInput(polygon)
        clean.Update()
        polygon = clean.GetOutput()
        
        pruned = self.prune(polygon)
        pruned = self.prune(pruned)
        polybound = chaste.interfaces.vtk_tools.vtk_tools.vtk_lines_to_polylines(pruned)
        
        smoother = vmtkscripts.vmtkCenterlineResampling()
        smoother.Centerlines = polybound
        smoother.length = 15.0
        smoother.Execute()
        self.network = smoother.Centerlines
        
        boundary_point_label = vtk.vtkFloatArray()
        boundary_point_label.SetNumberOfComponents(1)
        boundary_point_label.SetName("Radius")    
        
        locator = vtk.vtkKdTreePointLocator()
        locator.SetDataSet(self.surface) 
        vtk_numPoints = self.network.GetNumberOfPoints()    
        vtk_points = self.network.GetPoints() 
        surf_points = self.surface.GetPoints() 
        
        for idx in range(vtk_numPoints):
            loc = vtk_points.GetPoint(idx)
            closest = locator.FindClosestPoint(loc)
            surf_closest = surf_points.GetPoint(closest)
            radius = np.linalg.norm(np.array(loc) - np.array(surf_closest))/2.0
            boundary_point_label.InsertNextTupleValue((float(radius),))  
        
        self.network.GetPointData().AddArray(boundary_point_label)  

    def prune(self, surface):
        
        numPoints = surface.GetNumberOfPoints()    
        vtkpoints = surface.GetPoints()  
        connectivity = []
        points = [] 
        for i in range(numPoints):
            points.append(vtkpoints.GetPoint(i))
            connectivity.append([])
            
        numCells = surface.GetNumberOfLines()  
        cellArray = surface.GetLines()
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
            
        c1_points = np.zeros(numPoints)
        locator = vtk.vtkKdTreePointLocator()
        locator.SetDataSet(surface)
        near_seeds = np.zeros(numPoints)
        for idx, eachLoc in enumerate(self.start_points):
            probe_loc = np.array((eachLoc[0], eachLoc[1], 0.0))
            closest_id = locator.FindClosestPoint(probe_loc)
            near_seeds[closest_id] = 1.0
        
        for idx in range(points):
            num_segs = len(connectivity[idx])
            if num_segs == 1 and near_seeds[idx]==0:
                c1_points[idx] = 1
                current_point_id = idx
                previous_point_id = idx
                found_branch = False
                while not found_branch:
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
                        next_edge = edges[connectivity[idx][0]]
                              
                    if next_edge[0]== current_point_id:
                        next_point_id = next_edge[1]  
                    else:
                        next_point_id = next_edge[0]     
                              
                    if len(connectivity[next_point_id]) == 2:
                        c1_points[next_point_id] = 1
                        previous_point_id = current_point_id
                        current_point_id = next_point_id
                    else:
                        found_branch = True
                        break
            
        boundary_point_label = vtk.vtkFloatArray()
        boundary_point_label.SetNumberOfComponents(1)
        boundary_point_label.SetName("PointRemoveLabel")     
        
        for idx in range(c1_points):
            boundary_point_label.InsertNextTupleValue((float(c1_points[idx]),))  
            
        surface.GetPointData().AddArray(boundary_point_label)  
        
        threshold = vtk.vtkThreshold()
        threshold.SetInput(surface)
        threshold.SetInputArrayToProcess(0, 0, 0, "vtkDataObject::FIELD_ASSOCIATION_POINTS", "PointRemoveLabel")
        threshold.ThresholdBetween(0, 0)
        threshold.Update()
            
        polysurface = vtk.vtkGeometryFilter()
        polysurface.SetInputConnection(threshold.GetOutputPort())
        polysurface.Update()
        
        return polysurface.GetOutput()
        
    def get_output(self):
        return self.network
    
"""
Convert meshes between Tri/Tetgen and Dolfin formats
"""

from dolfin import *
import numpy as np
import chaste.population.vessel
    
class Extract2d():
    
    def __init__(self, mesh, boundary_markers, boundary_label):
        self.mesh = mesh
        self.boundary_markers = boundary_markers
        self.boundary_label = boundary_label
        
    def generate(self):
        
        # Get the centreline mesh
        facet_count = 0
        centre_facets = []
        unique_verts = []
        vert_coords = []
        
        for eachFacet in facets(self.mesh):
            if self.boundary_markers[facet_count] == self.boundary_label:
                v0 = Vertex(self.mesh, eachFacet.entities(0)[0])
                v1 = Vertex(self.mesh, eachFacet.entities(0)[1])
                centre_facets.append((v0, v1))
                
                if v0.index() not in unique_verts:
                    unique_verts.append(v0.index())
                    vert_coords.append(v0)
                if v1.index() not in unique_verts:
                    unique_verts.append(v1.index())
                    vert_coords.append(v1)   
            facet_count += 1      
            
        # construct the points 
        editor = MeshEditor()
        centre_mesh = Mesh()
        
        editor.open(centre_mesh, 1, 2)
        editor.init_vertices(len(unique_verts))
        editor.init_cells(len(centre_facets))
        
        nodes = []
        network = chaste.population.vessel.VascularNetwork()
        for idx, eachPoint in enumerate(vert_coords):
            editor.add_vertex(idx, np.array((eachPoint.x(0), eachPoint.x(1)), dtype=np.float))
            nodes.append(chaste.population.vessel.VascularNode(np.array((eachPoint.x(0), eachPoint.x(1), 0.0))))
        for idx, eachFacet in enumerate(centre_facets):
            id1 = unique_verts.index(eachFacet[0].index())
            id2 = unique_verts.index(eachFacet[1].index())
            editor.add_cell(idx, np.array((id1, id2), dtype=np.uintp))
            network.addVessel(chaste.population.vessel.Vessel([nodes[id1], nodes[id2]]))
            
        # Add extensions
        editor.close()
        centre_mesh.init()
        centre_mesh.order()
        
        return centre_mesh, network
    
class Extract3d():
    
    def __init__(self, mesh, edge_markers, edge_label):
        self.mesh = mesh
        self.edge_markers = edge_markers
        self.edge_label = edge_label
        
    def generate(self):
        
        # Get the centreline mesh
        edge_count = 0
        centre_edges = []
        unique_verts = []
        vert_coords = []
        
        for eachEdge in edges(self.mesh):
            if self.edge_markers[edge_count] == self.edge_label:
                v0 = Vertex(self.mesh, eachEdge.entities(0)[0])
                v1 = Vertex(self.mesh, eachEdge.entities(0)[1])
                centre_edges.append((v0, v1))
                
                if v0.index() not in unique_verts:
                    unique_verts.append(v0.index())
                    vert_coords.append(v0)
                if v1.index() not in unique_verts:
                    unique_verts.append(v1.index())
                    vert_coords.append(v1)   
            edge_count += 1      
            
        # construct the points 
        editor = MeshEditor()
        centre_mesh = Mesh()
        
        editor.open(centre_mesh, 1, 3)
        editor.init_vertices(len(unique_verts))
        editor.init_cells(len(centre_edges))
        
        nodes = []
        network = chaste.population.vessel.VascularNetwork()
        
        for idx, eachPoint in enumerate(vert_coords):
            editor.add_vertex(idx, np.array((eachPoint.x(0), eachPoint.x(1), eachPoint.x(2)), dtype=np.float))
            nodes.append(chaste.population.vessel.VascularNode(np.array((eachPoint.x(0), eachPoint.x(1), eachPoint.x(2)))))
        for idx, eachEdge in enumerate(centre_edges):
            id1 = unique_verts.index(eachEdge[0].index())
            id2 = unique_verts.index(eachEdge[1].index())
            editor.add_cell(idx, np.array((id1, id2), dtype=np.uintp))
            network.addVessel(chaste.population.vessel.Vessel([nodes[id1], nodes[id2]]))
            
        # Add extensions
        editor.close()
        centre_mesh.init()
        centre_mesh.order()
        
        return centre_mesh, network