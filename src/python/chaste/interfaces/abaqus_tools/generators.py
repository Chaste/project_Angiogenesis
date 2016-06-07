import os
import cPickle as pickle
import numpy as np
from abaqus import *
from abaqusConstants import *
import mesh
import step

def add_vessel(model, reference_part, index, start, end, radius):
    
    # Get the length
    length = np.linalg.norm(end - start)

    # Copy the reference part to a new part
    part_name = 'Vessel-' + str(index)
    p = model.Part(name=part_name, objectToCopy=reference_part)

#   # Edit the sketch
    s = p.features['Solid revolve-1'].sketch
    model.ConstrainedSketch(name='__edit__', objectToCopy=s)
    s = model.sketches['__edit__']
    s.parameters['radius'].setValues(expression=str(radius))
    s.parameters['length'].setValues(expression=str(length))
    p.features['Solid revolve-1'].setValues(sketch=s)
    p.regenerate()

    # Add the new part to the assembly, we need to translate and rotate it into the correct position
    midpoint = (start + end) / 2.0
    direction = (end - start) / length
    y_axis = np.array((0.0, 1.0, 0.0))
    rotation_angle = np.arccos(np.dot(y_axis, direction))
    rotation_axis = np.cross(y_axis, direction)
    
    a = model.rootAssembly
    a.regenerate()
    a.Instance(name=part_name, part=p, dependent=ON)
    a.rotate(instanceList=(part_name, ), axisPoint=(0.0, 0.0, 0.0), 
             axisDirection=tuple(rotation_axis), angle=180.0*(rotation_angle/np.pi))
    a.translate(instanceList=(part_name, ), vector=midpoint)
    
def modify_domain(bbox, model, reference_part):
    
        # Copy the reference part to a new part
    part_name = 'Domain-Temp'
    p = model.Part(name=part_name, objectToCopy=reference_part)
    
    #   # Edit the sketch
    s = p.features['Solid extrude-1'].sketch
    model.ConstrainedSketch(name='__edit__', objectToCopy=s)
    s = model.sketches['__edit__']
    s.parameters['ydim'].setValues(expression=str(bbox[3] - bbox[2]))
    s.parameters['xdim'].setValues(expression=str(bbox[1] - bbox[0]))
    p.features['Solid extrude-1'].setValues(sketch=s)
    p.regenerate()
    p.features['Solid extrude-1'].setValues(depth=bbox[5] - bbox[4])
    p.regenerate()
    
    a = model.rootAssembly
    a.regenerate()
    a.Instance(name=part_name, part=p, dependent=ON)
    a.translate(instanceList=(part_name, ), vector=(0.0, 0.0, -(bbox[5] - bbox[4])/2.0))
    
def create_material_and_sections(m, p, name = "Tissue"):
    m.Material(name=name)
    m.materials[name].Conductivity(table=((1.0, ), ))
    m.HomogeneousSolidSection(name=name, material=name)
    region = p.Set(cells=p.cells[0:1], name=name)
    p.SectionAssignment(region=region, sectionName=name)
    
def mesh_part(p, size=0.1):
    p.setMeshControls(regions=p.cells, elemShape=TET, technique=FREE)
    p.seedPart(size=size, deviationFactor=0.1, minSizeFactor=0.1)
    elemType1 = mesh.ElemType(elemCode=DC3D4, elemLibrary=STANDARD)
    p.setElementType(regions=(p.cells,), elemTypes=(elemType1,))
    p.generateMesh()

# Read in the vessel network and domain data
workdir = os.getcwd()
with open(workdir + "/model_info.pickle", 'rb') as handle:
    data  = pickle.load(handle)
network_data = data[0]
domain_data = data[1]
boundary_type = data[2]
simulation_name = data[3]

# Open the cae file containing the primitive geometry and make a copy of the model                         
openMdb(pathName=workdir + '/Model.cae')
mdb.Model(name=simulation_name, objectToCopy=mdb.models['Primitive'])
m = mdb.models[simulation_name]

# Create individual parts for each vessel in the network
p = m.parts['Vessel']
for idx,eachVessel in enumerate(network_data):
    start = eachVessel["Start Node Location"]
    end = eachVessel["End Node Location"]
    radius = eachVessel["Radius"]
    add_vessel(m, p, idx, np.array(start), np.array(end), radius)
    
# Create the domain
p = m.parts['Domain']
modify_domain(domain_data['BBox'], m, p)

# Do the merge and cut
a = m.rootAssembly
instances = []
for idx in range(len(network_data)):
    instances.append(a.instances['Vessel-' + str(idx)])

a.InstanceFromBooleanMerge(name='Merged', instances=instances, originalInstances=SUPPRESS, domain=GEOMETRY)
a.InstanceFromBooleanCut(name='SimDomain', instanceToBeCut=a.instances['Domain-Temp'], cuttingInstances=(a.instances['Merged-1'], ), 
        originalInstances=SUPPRESS)


# Set up material, section, mesh, step
p = m.parts['SimDomain']
create_material_and_sections(m, p)
mesh_part(p, (domain_data['BBox'][3]-domain_data['BBox'][2])*0.1)

# Set up simulation steps and output
m.HeatTransferStep(name='Step-1', previous='Initial', response=STEADY_STATE, amplitude=RAMP)
m.fieldOutputRequests['F-Output-1'].setValues(variables=('NT', 'TEMP', 'HFL', 'RFL', 'EVOL'))

# Add the body force
x_axis = np.array((1.0, 0.0, 0.0))
y_axis = np.array((0.0, 1.0, 0.0))
z_axis = np.array((0.0, 0.0, 1.0))

faces = a.instances['SimDomain-1'].faces
for idx, eachFace in enumerate(a.instances['SimDomain-1'].faces):
    not_cube = True
    face_normal = np.array(eachFace.getNormal())
    if np.array_equal(face_normal, x_axis):
        not_cube = False
    if np.array_equal(face_normal, -x_axis):
        not_cube = False    
    if np.array_equal(face_normal, y_axis):
        not_cube = False
    if np.array_equal(face_normal, -y_axis):
        not_cube = False
    if np.array_equal(face_normal, z_axis):
        not_cube = False
    if np.array_equal(face_normal, -z_axis):
        not_cube = False
    
    # Need to explicitly slice the Abaqus face array, rather than create a new list
    if not_cube:
        if idx==0:
            bound_faces = faces[eachFace.index:eachFace.index+1]
        else:
            bound_faces = bound_faces + faces[eachFace.index:eachFace.index+1]

a.Set(cells=a.instances['SimDomain-1'].cells[0:], name='BodyForceSet')
a.Surface(side1Faces=bound_faces, name='BcSurface')
a.Set(faces=bound_faces, name='BcSurfaceSet')
region = a.sets['BodyForceSet']
m.BodyHeatFlux(name='Load-1', createStepName='Step-1', region=region, magnitude=1.0, distributionType=USER_DEFINED)

if boundary_type == "Dirichlet":
    region = a.sets['BcSurfaceSet']
    m.TemperatureBC(name='BC-1', createStepName='Step-1', region=region, fixed=OFF, distributionType=UNIFORM, fieldName='', 
        magnitude=1.0, amplitude=UNSET)
else:
    region = a.surfaces['BcSurface']
    m.SurfaceHeatFlux(name='Load-2', createStepName='Step-1', region=region, magnitude=1.0, distributionType=USER_DEFINED)
 
# Create the job files
mdb.Job(name=simulation_name, model=simulation_name)
mdb.jobs[simulation_name].writeInput(consistencyChecking=OFF)
