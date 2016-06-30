from meshpy.triangle import MeshInfo, build
import chaste.utility.bases as bases
import chaste.mesh.converters as converters

class ChasteGeometryMesher2d(bases.SimpleIOBase):
    
    def __init__(self):
        super(ChasteGeometryMesher2d, self).__init__()
        
        self.mesh_size = 10.0
        self.region1_mesh_size = 10.0
        self.region2_mesh_size = 10.0
        self.region1_points = None
        self.region2_points = None
        self.holes = None
        
    def update(self):
        
        # convert to vtk
#        vtk_geometry = self.input.GetVtk(True)
        vtk_geometry = self.input
      
        # convert vtk to MeshPy Tri format
        converter = converters.VtkToTriMesh()
        converter.input = vtk_geometry
        converter.update()
        points, edges = converter.output
        
        print points, edges
        
        # Do the meshing with triangle
        mesh_info = MeshInfo()
        mesh_info.set_points(points)
        mesh_info.set_facets(edges)
        
        # Allow for two different regions and holes
        if self.region1_points is not None:
            total_regions = len(self.region1_points)
            if self.region2_points is not None:
                total_regions += len(self.region2_points)
            mesh_info.regions.resize(total_regions)
            for idx, eachPoint in enumerate(self.region1_points):
                mesh_info.regions[idx] = [eachPoint[0],eachPoint[1],1,self.region1_mesh_size]
            if self.region2_points is not None:
                for idx, eachPoint in enumerate(self.region2_points):
                    mesh_info.regions[idx + len(self.region1_points)] = [eachPoint[0],eachPoint[1],2,self.region2_mesh_size]
            if self.holes is not None:
                mesh_info.holes.resize(len(self.holes))
                for idx, eachHole in enumerate(self.holes):
                    mesh_info.holes[idx] = [eachHole[0], eachHole[1]]             
            data = build(mesh_info, volume_constraints=True, attributes=True, generate_faces=True)
            self.output = [data.points, data.elements]
        else:
            if self.holes is not None:
                mesh_info.holes.resize(len(self.holes))
                for idx, eachHole in enumerate(self.holes):
                    mesh_info.holes[idx] = [eachHole[0], eachHole[1]]
            data = build(mesh_info, max_volume = self.region1_mesh_size)
            self.output = [data.points, data.elements]
        return self.output
        
class EmbeddedChasteGeometryMesher2d(bases.SimpleIOBase):
    
    def __init__(self):
        super(EmbeddedChasteGeometryMesher2d, self).__init__()
        self.mesh_size = 10.0
        self.region1_mesh_size = 10.0
        self.region2_mesh_size = 10.0
        self.region1_points = None
        self.region2_points = None
        self.holes = None
        self.embedded_geometry = None
    
    def update(self):    
        
        # convert geometries to vtk
        outer_vtk_geometry = self.input.GetVtk(True)
        inner_vtk_geometry = self.embedded_geometry.GetVtk(True)
        
        # convert vtk to MeshPy Tri format
        converter = converters.VtkToTriMesh()
        converter.input = outer_vtk_geometry
        converter.update()
        points, edges = converter.output
        
        converter.input = inner_vtk_geometry
        converter.update()
        inner_points, inner_edges = converter.output   
        inner_edges = [[x[0]+len(points), x[1]+len(points)] for x in inner_edges]     
        points.extend(inner_points)
        edges.extend(inner_edges)
            
        # Do the meshing with triangle
        mesh_info = MeshInfo()
        mesh_info.set_points(points)
        mesh_info.set_facets(edges)
        
        if self.region1_points is not None:
            total_regions = len(self.region1_points)
            if self.region2_points is not None:
                total_regions += len(self.region2_points)
            
            mesh_info.regions.resize(total_regions)
            for idx, eachPoint in enumerate(self.region1_points):
                mesh_info.regions[idx] = [eachPoint[0],eachPoint[1],1,self.region1_mesh_size]
            if self.region2_points is not None:
                for idx, eachPoint in enumerate(self.region2_points):
                    mesh_info.regions[idx + len(self.region1_points)] = [eachPoint[0],eachPoint[1],2,self.region2_mesh_size]
            if self.holes is not None:
                mesh_info.holes.resize(len(self.holes))
                for idx, eachHole in enumerate(self.holes):
                    mesh_info.holes[idx] = [eachHole[0], eachHole[1]]       
            self.output = build(mesh_info, volume_constraints=True, attributes=True, generate_faces=True)
        else:
            
            if self.holes is not None:
                mesh_info.holes.resize(len(self.holes))
                for idx, eachHole in enumerate(self.holes):
                    mesh_info.holes[idx] = [eachHole[0], eachHole[1]]
            self.output = build(mesh_info, max_volume = self.region1_mesh_size)
        return self.output
    
import collections
import vtk
import chaste.mesh
#import code.geometry.converters
import subprocess
import pandas
try:
   import cPickle as pickle
except:
   import pickle
    
class Mesher3d():
    
    def __init__(self, work_dir, domain = None, network = None, region_markers = None, surface = True, mesh_size= 100.0, vtk_surface = None):
        self.domain = domain
        self.network = network
        self.region_markers = region_markers
        self.surface = surface
        self.mesh_size = mesh_size
        self.vtk_surface = vtk_surface
        self.work_dir = work_dir
        
    def generate(self):
        
        if not self.surface:
            self.domain.AddVesselNetwork(self.network, False)
            mesh = chaste.mesh.HybridMesh()
            mesh.GenerateFromPart(self.domain, self.mesh_size)
            return mesh
        
        elif self.surface and not self.vtk_surface:
            
            # generate the stl
            #net_converter = code.geometry.converters.NetworkToGeometry(self.network)
            #network_geometry = net_converter.generate()
            
            #converter = code.geometry.converters.PartToGeometry(self.domain)
            #domain_geometry = self.domain
            #full_geometry = domain_geometry.fuse(network_geometry)
            #full_geometry.exportStl(self.work_dir + "/merge.stl")
            
            mesh_command = "cd " + self.work_dir + "; tetgen -pq1.5kg merge.stl"
            p =  subprocess.Popen(mesh_command, shell=True)
            p.wait()
              
            # Add the regions to the smesh file
            if self.region_markers is not None:
                readFile = open(self.work_dir + "/merge.1.smesh")
                lines = readFile.readlines()
                readFile.close()
                  
                total_regions = 0
                for eachRegion in self.region_markers:
                    if isinstance(eachRegion[0], collections.Sequence):
                        total_regions += len(eachRegion[0])
                    else:
                        total_regions += 1
                  
                w = open(self.work_dir + "/merge.1.smesh",'w')
                w.writelines([item for item in lines[:-2]])
                w.write(str(total_regions) + " \n")
                
                counter = 0
                for eachRegion in self.region_markers:
                    for idx, eachPoint in enumerate(eachRegion[0]):
                        w.write(str(counter) + " " + str(eachPoint[0]) + " " + str(eachPoint[1]) + " " + str(eachPoint[2]) + " "
                                 + str(eachRegion[1]) + " " + str(eachRegion[2]) +  "\n")
                        counter += 1
                        
                w.write(" \n")
                w.close()
                
                mesh_command = "tetgen -pkgAq1.5 " + self.work_dir + "/merge.1.node " + self.work_dir + "/merge.1.smesh"
                p =  subprocess.Popen(mesh_command, shell=True)
                p.wait()
                
                convert_command = "dolfin-convert " + self.work_dir + "/merge.2.mesh " + self.work_dir + "/dolfin_mesh.xml"
                print subprocess.Popen(convert_command, shell=True, stdout=subprocess.PIPE).stdout.read()
                elem_data = pandas.read_csv(self.work_dir + "/merge.2.ele", delim_whitespace=True, header=None, error_bad_lines=False, skiprows=1)
                
                with open(self.work_dir + "/physical_labels.pickle", 'wb') as handle:
                    pickle.dump(elem_data[5], handle) 
            else:
                convert_command = "dolfin-convert " + self.work_dir + "/merge.1.mesh " + self.work_dir + "/dolfin_mesh.xml"
                print subprocess.Popen(convert_command, shell=True, stdout=subprocess.PIPE).stdout.read()