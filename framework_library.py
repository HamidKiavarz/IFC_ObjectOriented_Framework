import numpy as np
import ifcopenshell
from shapely.geometry import Point
from shapely.geometry import LineString
from operator import itemgetter
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math
import OCC.Core.gp
import OCC.Core.Geom
from OCC.Core  import GProp, BRepGProp
import OCC.Core.Units
import OCC.Core.Bnd
import OCC.Core.BRepBndLib
import OCC.Core.BRep
import OCC.Core.BRepPrimAPI
import OCC.Core.BRepAlgoAPI
import OCC.Core.BRepBuilderAPI
import OCC.Core.GProp
import OCC.Core.BRepGProp
import ifcopenshell
import ifcopenshell.geom
from OCC.Core.TopoDS import (topods, TopoDS_Wire, TopoDS_Vertex, TopoDS_Edge,
                             TopoDS_Face, TopoDS_Shell, TopoDS_Solid, TopoDS_Shape,
                             TopoDS_Compound, TopoDS_CompSolid, topods_Edge,
                             topods_Vertex, TopoDS_Iterator)
from IFC.core_topology_traverse import Topo
from OCC.Core.BRep import BRep_Tool_Surface
from OCC.Core.BRep import BRep_Tool
from OCC.Display.SimpleGui import init_display
import IFC.Entities as model
import IFC.config as config
def norm2(X):
    return np.sqrt(np.sum(X ** 2))
def cmp_floats(a,b, atol=1e-12):
    return abs(a-b) < atol
def magnitude(vector):
    return np.sqrt(np.dot(np.array(vector), np.array(vector)))
def norm(vector):
    return np.array(vector)/magnitude(np.array(vector))
def normalized(X):
	return X / norm2(X)

def get_plane_plane_intersection_line(p1, p2):
    u = normalized(np.cross(p1[:-1], p2[:-1]))
    m = np.array((p1[:-1], p2[:-1], u))
    x = np.array((-p1[-1], -p2[-1], 0.))
    return u, np.linalg.solve(m, x)

def is_points_on_face(frmw_face, points):
    xmin_t = float("{:3f}".format(frmw_face.point_min.X))
    ymin_t = float("{:3f}".format(frmw_face.point_min.Y))
    zmin_t = float("{:3f}".format(frmw_face.point_min.Z))
    xmax_t = float("{:3f}".format(frmw_face.point_max.X))
    ymax_t = float("{:3f}".format(frmw_face.point_max.Y))
    zmax_t = float("{:3f}".format(frmw_face.point_max.Z))
    on_face = 0
    for p in points:
        x = float("{:3f}".format(p.X))
        y = float("{:3f}".format(p.Y))
        z = float("{:3f}".format(p.Z))
        if (xmin_t-0.5 <= x and x <= xmax_t+0.5):
            if (ymin_t-0.5 <= y and y <= ymax_t+0.5):
                # if(zmin_t <= z and z <= zmax_t):
                on_face += 1
    if on_face == 2:
        return True
    else:
        return False
def polygon_area(ptn1,ptn2,h):
    p1 = np.array(ptn1)
    p2 = np.array(ptn2)
    squared_dist = abs(np.sum((p1[:-1]-p2[:-1])**2, axis=0))
    dist = np.sqrt(squared_dist)
    area = dist * h
    return float("{:3f}".format(area))
def make_plane(input_face):
    face_bbox = OCC.Core.Bnd.Bnd_Box()
    OCC.Core.BRepBndLib.brepbndlib_Add(input_face,face_bbox) 
    face_center = ifcopenshell.geom.utils.get_bounding_box_center(face_bbox).XYZ()
    surf = BRep_Tool_Surface(input_face)
    pln = OCC.Core.Geom.Handle_Geom_Plane_DownCast(surf)
    normal = pln.Axis().Direction().XYZ()
    d = -np.sum(face_center*normal)
    return [float("{:3f}".format(normal.X())), float("{:3f}".format(normal.Y())), float("{:3f}".format(normal.Z())), float("{:3f}".format(d))]
def plane_info(face_center,face_normal):
    point  = np.array([face_center.X(),face_center.Y(),face_center.Z()])
    normal = np.array([float("{:3f}".format(face_normal.X())),float("{:3f}".format(face_normal.Y())),float("{:3f}".format(face_normal.Z()))])
    # a plane is a*x+b*y+c*z+d=0
    # [a,b,c] is the normal. Thus, we have to calculate
    d = -np.sum(point*normal)# dot product
    eq = '%2.1fX + %2.1fY + %2.1fZ  = %2.2f' %(normal[0],normal[1],normal[2],d*-1)
    return [eq,normal,float("{:3f}".format(d))]  
def space_bbx_center(space,settings):
    space_bbox = OCC.Core.Bnd.Bnd_Box()
    space_shape = ifcopenshell.geom.create_shape(settings, space).geometry
    OCC.Core.BRepBndLib.brepbndlib_Add(space_shape,space_bbox) 
    space_bbox_center = ifcopenshell.geom.utils.get_bounding_box_center(space_bbox)
    return space_bbox_center
def distance_from_plane(point, plane_param):  
    x = point[0]
    y = point[1]
    z = point[2]
    a = plane_param[0]
    b = plane_param[1]
    c = plane_param[2]
    d = plane_param[3]

    d = abs((a * x + b * y + c * z + d))  
    e = (math.sqrt(a * a + b * b + c * c)) 
    return d/e
def project_point_on_plane(q, p, n):
    ### this method calculate the projected of given point on given plane
    pq = q - p
    dot = pq.Dot(n)
    dotn  = np.array([n.X(),n.Y(),n.Z()]) * dot
    q_proj = np.array([q.X(),q.Y(),q.Z()]) - dotn
    return q_proj
def calculate_face_area(face):
    props = GProp.GProp_GProps()
    BRepGProp.brepgprop_SurfaceProperties(face, props)
    surface_area = props.Mass()
    return float("{:3f}".format(surface_area))
def normalize(x):
    x = np.asarray(x)
    return (x - x.min()) / (np.ptp(x))
def PolygonArea(corners):
    n = len(corners) # of corners
    area = 0.0
    for i in range(n):
        j = (i + 1) % n
        area += corners[i][0] * corners[j][1]
        area -= corners[j][0] * corners[i][1]
    area = abs(area) / 2.0
    return area
def find_face_toward_center(space_bbx_center,face_obj_list):
    face_toward_center_info = []
    face_final_list = []
    
    face_obj_sorted_list = sorted(face_obj_list, key=lambda x:x[1])
    ##  select two main spaces (the larg face area of wall)
    face_final_list.append(face_obj_sorted_list[len(face_obj_sorted_list) - 2])
    face_final_list.append(face_obj_sorted_list[len(face_obj_sorted_list) - 1])
    distance = 1000
    dist = 0
    for item in face_final_list:
        face = item[0]
        surf = BRep_Tool_Surface(face)
        plane = OCC.Core.Geom.Handle_Geom_Plane_DownCast(surf)
        face_normal = plane.Axis().Direction().XYZ()
        # face_normal_vector = np.array([face_normal.X(),face_normal.Y(),face_normal.Z()])
        face_normal_vector = np.array([float("{:3f}".format(face_normal.X())),float("{:3f}".format(face_normal.Y())),float("{:3f}".format(face_normal.Z()))])
        face_bbox = OCC.Core.Bnd.Bnd_Box()
        OCC.Core.BRepBndLib.brepbndlib_Add(face,face_bbox) 
        space_center = np.array([space_bbx_center.X(),space_bbx_center.Y(),space_bbx_center.Z()])
        face_center = ifcopenshell.geom.utils.get_bounding_box_center(face_bbox).XYZ()
        prj_space_point = project_point_on_plane(space_bbx_center,face_center,face_normal)
        wall_space_vector = space_center - prj_space_point
        vector1 = wall_space_vector /np.linalg.norm(wall_space_vector)
        vector2 = face_normal_vector /np.linalg.norm(face_normal_vector)
        dot_product = np.dot(vector1, vector2)
        plane_d = plane_info(face_center,face_normal)[2]
        dist = distance_from_plane (space_center,[face_normal_vector[0],face_normal_vector[1],face_normal_vector[2],plane_d])
        if dist < distance:
            distance = dist
            face_toward_center_info = [face,item[1],face_normal_vector]
    return face_toward_center_info
def shape_display(shape):
    display, start_display, add_menu, add_function_to_menu = init_display()
    display.DisplayShape(shape, update=True)
    # display.SetSelectionModeEdge() # switch to edge selection mode
    display.SetSelectionModeFace()
    display.register_select_callback(line_clicked)

    start_display()
def line_clicked(shp, *kwargs):
    """ This function is called whenever a line is selected"""
    for shape in shp: # this should be a TopoDS_Edge
        props = GProp.GProp_GProps()
        BRepGProp.brepgprop_SurfaceProperties(shape, props)
        # BRepGProp.brepgprop_LinearProperties(target_face,l_props)
        face_area = props.Mass()
        print("shape selected: ", face_area)
def extract_internal_face(ifc_space,wall):
    wall_shape = None
    face_obj = model.Face()
    face_obj.ifc_face = ifc_space
    face_obj.name = ifc_space.Name + '_' + wall.Tag
    settings = ifcopenshell.geom.settings()
    settings.set(settings.USE_PYTHON_OPENCASCADE, True)
    wall_bbox = OCC.Core.Bnd.Bnd_Box()
    
    # Initialize a graphical display window
    try:
        wall_shape = ifcopenshell.geom.create_shape(settings, wall).geometry 
    except:
        wall_shape = None
        face_obj = None
        # print("Cannot create wall shape %s" %wall.GlobalId)
    # from OCC.Display.SimpleGui import init_display
    # shape_display(wall_shape)

    #add shape to bbox . It could be several shapes (walls of a room and space)
    # OCC.Core.BRepBndLib.brepbndlib_Add(wall_shape,wall_bbox) 
    
    ### calculate bbox size
    # xmin, ymin, zmin, xmax, ymax, zmax = wall_bbox.Get()
    # dx = xmax - xmin
    # dy = ymax - ymin
    # dz = zmax - zmin

    # Calculate the center/average of the bounding box of wall and space
    # wall_bbox_center = ifcopenshell.geom.utils.get_bounding_box_center(wall_bbox)
    
    ### Extracting faces
    if wall_shape != None:
        faces = list(Topo(wall_shape).faces())
        face_count = Topo(wall_shape).number_of_faces()
        space_bbx_cntr = space_bbx_center(ifc_space,settings).XYZ()
        
        face_obj_list = []
        for face in faces:
            face_obj_list.append([face,calculate_face_area(face)])
        face_toward_cntr_info = find_face_toward_center(space_bbx_cntr,face_obj_list)
        # face toward center(ftc)
        ftc = face_toward_cntr_info[0]        
        
        face_bbox = OCC.Core.Bnd.Bnd_Box()
        OCC.Core.BRepBndLib.brepbndlib_Add(ftc,face_bbox) 
        # xmin, ymin, zmin, xmax, ymax, zmax = face_bbox.Get()
        min_max_coords = face_bbox.Get()
        face_obj.wall = wall
        face_obj.point_min = model.Point(float("{:3f}".format(face_bbox.Get()[0])),float("{:3f}".format(face_bbox.Get()[1])), float("{:3f}".format(face_bbox.Get()[2])))
        face_obj.point_max = model.Point(float("{:3f}".format(face_bbox.Get()[3])),float("{:3f}".format(face_bbox.Get()[4])), float("{:3f}".format(face_bbox.Get()[5])))
        face_obj.normal_vector = face_toward_cntr_info[2]
        face_obj.ifc_face = ftc
        ### wal.guid,face object, min-max of face, face normal
        # return [wall.GlobalId,ftc,min_max_coords,face_toward_cntr_info[2]]
    return face_obj
def extract_space_minmax(ifc_space):
    settings = ifcopenshell.geom.settings()
    settings.set(settings.USE_PYTHON_OPENCASCADE, True)
    space_bbox = OCC.Core.Bnd.Bnd_Box()
    space_shape = ifcopenshell.geom.create_shape(settings, ifc_space).geometry 
    OCC.Core.BRepBndLib.brepbndlib_Add(space_shape,space_bbox) 
    xmin, ymin, zmin, xmax, ymax, zmax = space_bbox.Get()
    return([[xmin, ymin, zmin],[xmax, ymax, zmax]])
def find_floor_wall_faces(floor_number,all_ifc_spaces):
    # ifc_file = ifcopenshell.open(config.path)
    # all_ifc_spaces = ifc_file.by_type("IfcSpace")
    fmw_space_lst = []
    floor_wall_faces_dic = {}
    for ifc_s in all_ifc_spaces: # all spaces in the second floor
        if (ifc_s.Decomposes[0][4][2]).lower() == 'floor ' + str(floor_number):
            fmw_space_obj = model.Space(ifc_s)
            fmw_space_lst.append(fmw_space_obj)
    for fmw_s in fmw_space_lst: # create the dictionary of wall(key) and its faces(values)
            for fmw_w in fmw_s.fmw_walls:
                if fmw_w.guid not in floor_wall_faces_dic.keys():
                    floor_wall_faces_dic[fmw_w.guid] = []
                    floor_wall_faces_dic[fmw_w.guid].append(fmw_w.fmw_face)
                else:
                    floor_wall_faces_dic[fmw_w.guid].append(fmw_w.fmw_face)
    return floor_wall_faces_dic
def project_point_to_line(point, line):
    ''' project  2D point on 2D line 
    example: point = Point(0.2, 0.5)
    line = LineString([(0, 1), (1, 1)]) '''
    
    x = np.array(point.coords[0])

    u = np.array(line.coords[0])
    v = np.array(line.coords[len(line.coords)-1])

    n = v - u
    n /= np.linalg.norm(n, 2)

    P = u + n*np.dot(x - u, n)

    return P
   