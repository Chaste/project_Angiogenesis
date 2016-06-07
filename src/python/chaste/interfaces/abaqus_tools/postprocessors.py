from odbAccess import *
from abaqusConstants import *
odb = openOdb(path='result.odb')
lastFrame = odb.steps['Step-1'].frames[-1]

# Set up a volume dict
volumes = lastFrame.fieldOutputs['EVOL']
volumeValues = volumes.values
volumes_dict = {}
for v in volumeValues:
    volumes_dict[v.elementLabel] = v.data

pressures = lastFrame.fieldOutputs['TEMP'].getSubset(position=CENTROID)
pressureValues = pressures.values

f = open('processed_result.dat','w')
for idx, v in enumerate(pressureValues):
    f.write(str(volumes_dict[v.elementLabel]) + ", " + str(v.data) +'\n')
f.close()