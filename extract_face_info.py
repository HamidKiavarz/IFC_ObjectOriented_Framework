import numpy as np
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
def extract_internal_face(space,wall):
    face_obj = model.Face()
    
    settings = ifcopenshell.geom.settings()
    settings.set(settings.USE_PYTHON_OPENCASCADE, True)
    wall_bbox = OCC.Core.Bnd.Bnd_Box()
    
    # Initialize a graphical display window
    wall_shape = ifcopenshell.geom.create_shape(settings, wall).geometry 
    # from OCC.Display.SimpleGui import init_display
    # shape_display(wall_shape)

    #add shape to bbox . It could be several shapes (walls of a room and space)
    OCC.Core.BRepBndLib.brepbndlib_Add(wall_shape,wall_bbox) 
    
    ### calculate bbox size
    xmin, ymin, zmin, xmax, ymax, zmax = wall_bbox.Get()
    dx = xmax - xmin
    dy = ymax - ymin
    dz = zmax - zmin

    # Calculate the center/average of the bounding box of wall and space
    wall_bbox_center = ifcopenshell.geom.utils.get_bounding_box_center(wall_bbox)
    
    ### Extracting faces
    faces = list(Topo(wall_shape).faces())
    face_count = Topo(wall_shape).number_of_faces()
    space_bbx_cntr = space_bbx_center(space,settings).XYZ()
    
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






    # space_B203 = ifc_file['0BTBFw6f90Nfh9rP1dl_39']
    # space_A205 = ifc_file['2gRXFgjRn2HPE$YoDLX3FV']
    # space_A201 = ifc_file['0BTBFw6f90Nfh9rP1dlXri']
    # space_A204 = ifc_file['0BTBFw6f90Nfh9rP1dlXre']
    # wall = ifc_file['2O2Fr$t4X7Zf8NOew3FKau']
    # wall = ifc_file['2O2Fr$t4X7Zf8NOew3FLEJ']
    # wall = ifc_file['2O2Fr$t4X7Zf8NOew3FL9r']
    # wall = ifc_file['2O2Fr$t4X7Zf8NOew3FLPP']
    # wall = ifc_file['2O2Fr$t4X7Zf8NOew3FLMr']
    # wall = ifc_file['2O2Fr$t4X7Zf8NOew3FLQD']
    # wall = ifc_file['2O2Fr$t4X7Zf8NOew3FLR9']
    # wall = ifc_file['2O2Fr$t4X7Zf8NOew3FLOH']

    # print ("Bounding box center: %.2f %.2f %.2f" % (
    #     wall_bbox_center.X(), 
    #     wall_bbox_center.Y(),
    #     wall_bbox_center.Z()))
    ### ---------------------------------------------

    ''' https://www.mail-archive.com/search?l=pythonocc-users@gna.org&q=subject:%22%5C%5BPythonocc%5C-users%5C%5D+retrieve+vertices+of+an+edge%22&o=newest&f=1'''
    # vertics = list(Topo(target_face).vertices())
    # coords = []
    # for vertex  in vertics:
    #     pnt = BRep_Tool().Pnt(vertex)
        
    #     # print('X: %.3f Y: %.3f Z: %.3f' % (pnt.X(),pnt.Y(),pnt.Z()))
    #     coords.append([pnt.X(),pnt.Y(),pnt.Z()])
    
    # # Inverse normal
    
    # print(pln.Axis().Direction().XYZ())
    # ### point on plane
    # print(pln.Location().X())
    # print(pln.Location().Y())
    # print(pln.Location().Z())