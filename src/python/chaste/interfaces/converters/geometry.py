import chaste.utility.standalone_bases as bases

class ChastePartToCadShell(bases.SimpleIOBase):
    
    def __init__(self):
        super(ChastePartToCadShell, self).__init__()
        
        self.BaseModule = __import__("FreeCAD")
        self.PartModule = __import__("Part")
        
    def update(self):
        polygons = self.input.GetPolygons()
        faces = []
        for eachPolygon in polygons:
            locs = []
            verts = eachPolygon.GetVertices()
            for eachVert in verts:
                locs.append(tuple(eachVert.GetLocation()))
            locs.append(verts[0].GetLocation())
            poly = self.PartModule.makePolygon(locs)
            faces.append(self.PartModule.Face(poly))
            
        self.output = self.PartModule.makeShell(faces)
        return self.output