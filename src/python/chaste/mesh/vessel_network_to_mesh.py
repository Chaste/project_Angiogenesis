import vtk
import dolfin as df
import numpy as np
import chaste.utility.bases as bases
import chaste.mesh.converters as converters
import chaste.population.vessel
import chaste.image.network_to_image

class VesselNetworkMesher1d(bases.SimpleIOBase):
    
    def __init__(self):
        super(VesselNetworkMesher1d, self).__init__()
        
        self.mesh_size = 10.0
        self.su_subdivisions = 10
        
    def update(self):
        
        # Get the vtk representation
        writer = chaste.population.vessel.VtkVesselNetworkWriter()
        writer.SetVesselNetwork(self.input)
        vtk_rep = writer.GetOutput()
        
        # subdivide
        spline = vtk.vtkSplineFilter();
        spline.SetInput(vtk_rep)
#        spline.SetLength(self.mesh_size*10.0)
        spline.SetNumberOfSubdivisions(10)
        
        spline.Update();
        
        triangle = vtk.vtkTriangleFilter()
        triangle.SetInput(spline.GetOutput())
        triangle.Update()
        
        clean = vtk.vtkCleanPolyData()
        clean.SetInput(triangle.GetOutput())
        clean.Update()
        
        vtk_rep = clean.GetOutput();
        
        # convert vtk to MeshPy Tri format
        converter = converters.VtkToTriMesh()
        converter.input = vtk_rep
        converter.update()
        points, edges = converter.output
        
        editor = df.MeshEditor()
        centre_mesh = df.Mesh()
        
        editor.open(centre_mesh, 1, 2)
        editor.init_vertices(len(points))
        editor.init_cells(len(edges))

        for idx, eachPoint in enumerate(points):
            editor.add_vertex(idx, np.array((eachPoint[0], eachPoint[1]), dtype=np.float))
            
        for idx, eachCell in enumerate(edges):
            editor.add_cell(idx, np.array((eachCell[0], eachCell[1]), dtype=np.uintp))

        editor.close()
        centre_mesh.init()
        centre_mesh.order()
        
        file = df.File("/home/grogan/ABME16_Work/test/TestSingleVesselMeshing/mesh_1d.xml")
        file << centre_mesh
        return centre_mesh
    
class VesselNetworkMesher2d(bases.SimpleIOBase):
     
    def __init__(self):
        super(VesselNetworkMesher1d, self).__init__()
         
        self.mesh_size = 10.0
        self.su_subdivisions = 10
         
    def update(self):
         
        # Get the image representation, get the edges, smooth, mesh 2d
        image = chaste.image.network_to_image.network_to_image(self.input, workdir, filename, grid, radius)
         

        