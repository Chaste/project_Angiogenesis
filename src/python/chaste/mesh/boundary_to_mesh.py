from meshpy.triangle import MeshInfo, build
import chaste.utility.bases as bases
import chaste.mesh.converters as converters

class BoundaryToMesh2d(bases.SimpleIOBase):
    
    def __init__(self):
        super(BoundaryToMesh2d, self).__init__()
        self.mesh_size = 10.0
        self.holes = None
        
    def update(self):
        
        # convert vtk to MeshPy Tri format
        converter = converters.VtkToTriMesh()
        converter.input = self.input
        converter.update()
        points, edges = converter.output
        
        # Do the meshing with triangle
        mesh_info = MeshInfo()
        mesh_info.set_points(points)
        mesh_info.set_facets(edges)
        
        if self.holes is not None:
            mesh_info.holes.resize(len(self.holes))
            for idx, eachHole in enumerate(self.holes):
                mesh_info.holes[idx] = [eachHole[0], eachHole[1]]             
        self.output = build(mesh_info, volume_constraints=True, attributes=True, generate_faces=True)
            
        return self.output
    
