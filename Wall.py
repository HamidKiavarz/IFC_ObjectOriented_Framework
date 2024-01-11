import ifcopenshell
import ifcopenshell.util
import ifcopenshell.util.element
from ifcopenshell import geom
from OCC.Core.TopoDS import (topods, TopoDS_Wire, TopoDS_Edge,
                             TopoDS_Face, TopoDS_Shell, TopoDS_Solid, TopoDS_Shape,
                             TopoDS_Compound, TopoDS_CompSolid, topods_Edge,
                             topods_Vertex, TopoDS_Iterator)
from OCC.Core.BRep import BRep_Tool
import numpy as np
import IFC.framework_library as lib 
def extract_vertics(wall):
    points = []
    line  = wall.Representation.Representations[0].Items[0]
    for item in wall.Representation.Representations:
        for t in item:
            print(t)
    # for p in line:
    #     for cPoint in p:
    #         print(cPoint.Coordinates)

def read_geom_as_mesh(entity):
    settings = geom.settings()
    settings.set(settings.USE_PYTHON_OPENCASCADE, True)
    shape = ifcopenshell.geom.create_shape(settings, entity)
    # ios stands for IfcOpenShell
    ios_vertices = shape.geometry.verts
    ios_edges = shape.geometry.edges
    ios_faces = shape.geometry.faces
    points = list(ios_vertices[x:x + 3]  for x in range(0, len(ios_vertices), 3)) 

    vertices = np.array(points)
    x=vertices[:,0]
    y=vertices[:,1]
    z=vertices[:,2]

    edges = [ios_edges[i : i + 2] for i in range(0, len(ios_edges), 2)]
    faces = [tuple(ios_faces[i : i + 3]) for i in range(0, len(ios_faces), 3)]
    print ("This %s has been defined by %s vertices, %s edges and %s faces" % (entity.is_a(),len(vertices),len(edges), len(faces)))
    print(vertices)
if __name__ == "__main__":
    
    path = r'C:\Users\Malavan-PC\OneDrive - York University\Data\IFC\DuplexModel-IFC\Duplex_A_20110505.ifc'
    ifc_file = ifcopenshell.open(path)
    ifc_wall = ifc_file['2O2Fr$t4X7Zf8NOew3FLR9']
    all_walls_second_floors = ifc_file.by_type("IfcWall")
    ifc_walls_f2 =  []
    for w in all_walls_second_floors:
        if w.Decomposes[0][4][2] == 'Level ' + str(2):
            ifc_walls_f2.append(w)
    floor_wall_dics = lib.find_floor_wall_faces(2)
    # print(ifcopenshell.util.element.get_psets(wall))
    '''Some attributes are special, called "inverse attributes". 
    They happen when another element is referencing our element. 
    They can reference it for many reasons, like to define a relationship, 
    such as if they create a void in our wall, join our wall, or define a 
    quantity take-off value for our wall, among others. Just treat them like regular attributes'''
    # print(wall.IsDefinedBy)
    print(ifc_wall.attribute_type)
    wall_property_set = None
    # for definition in ifc_wall.IsDefinedBy:
    #     # To support IFC2X3, we need to filter our results.
    #     if definition.is_a('IfcRelDefinesByProperties'):
    #         property_set = definition.RelatingPropertyDefinition
    #         if(property_set.Name == 'Pset_WallCommon'): # Might return Pset_WallCommon
    #             wall_property_set = definition.RelatingPropertyDefinition
    #             break
    # for property in wall_property_set.HasProperties:
    #     if property.is_a('IfcPropertySingleValue') and property.Name == 'IsExternal':
    #         print(property.NominalValue.wrappedValue)
    # # read_geom_as_mesh(wall)
    # spaces = []
    # print(len(list(wall.ProvidesBoundaries)))
    #get all related spaces
    # for pb in wall.ProvidesBoundaries: # rbe --> ifcRelSpaceBoundary
    #     if not pb.RelatingSpace in spaces:
    #         spaces.append(pb.RelatingSpace)
    #         print(pb.RelatingSpace.Name)
    # print(len(spaces))





    