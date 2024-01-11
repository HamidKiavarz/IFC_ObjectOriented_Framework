import numpy as np
import ifcopenshell
import math
from shapely.geometry import Point as SHPoint
from shapely.geometry import LineString
import IFC.framework_library as lib
import IFC.config as config
class Space(object):
    guid = None
    tag = None
    name = None
    long_name = None
    is_external = None
    fmw_walls = None
    neighbour_ifc_spaces = None
    neighbour_fmw_spaces = None
    def __init__(self, ifc_space):
        self.ifc_space = ifc_space
        self.name = ifc_space.Name
        self.guid = ifc_space.GlobalId
        self.fmw_walls = []
        self.extract_topological_walls()
        self.neighbour_ifc_spaces = []
        self.neighbour_fmw_spaces = {}
        self.find_ifc_neighbour_spaces()
        # neighbour_fmw_spaces = self.find_fmw_neighbour_spaces(2)
    def get_by_guid(self,guid):
        path = r'C:\Users\Malavan-PC\OneDrive - York University\Data\IFC\DuplexModel-IFC\Duplex_A_20110505.ifc'
        ifc_file = ifcopenshell.open(path)
        space_obj = Space()
        space_obj.guid = Space.guid
        return space_obj
    def extract_topological_walls(self):
        space_ifc_walls =[]
        wall_tags = []
        for rbe in self.ifc_space.BoundedBy: # rbe --> ifcRelSpaceBoundary . Find space walls
            related_element = rbe.RelatedBuildingElement
            if related_element != None:
                if(related_element.is_a('IfcWall') or related_element.is_a('IfcCurtainWall')):
                    if not related_element.Tag in wall_tags:
                        wall_tags.append(related_element.Tag)
                        space_ifc_walls.append(related_element)
                    # cbp = rbe.Con`nectionGeometry.SurfaceOnRelatingElement
        for target_ifc_wall in space_ifc_walls:
            wall_frw_obj = Wall(target_ifc_wall)
            target_fmw_face = lib.extract_internal_face(self.ifc_space,target_ifc_wall)
            wall_frw_obj.fmw_face = target_fmw_face
            if wall_frw_obj.fmw_face != None:
                target_frmw_pln = lib.make_plane(target_fmw_face.ifc_face)
                intersected_points = []
                intersect_point_1 = None
                intersect_point_2 = None
                for w in space_ifc_walls:
                    if w != target_ifc_wall:
                        other_frm_face = lib.extract_internal_face(self.ifc_space,w)
                        if other_frm_face != None:
                            dot = np.dot(other_frm_face.normal_vector,target_fmw_face.normal_vector)
                            if dot == 0: ### perpendicular faces
                                other_frmw_pln = lib.make_plane(other_frm_face.ifc_face)
                                intersect_line = lib.get_plane_plane_intersection_line(target_frmw_pln,other_frmw_pln)
                                intersect_point_1 = Point(intersect_line[1][0],intersect_line[1][1],target_fmw_face.point_min.Z)
                                intersect_point_2 = Point(intersect_line[1][0],intersect_line[1][1],target_fmw_face.point_max.Z)

                                if (lib.is_points_on_face(other_frm_face, [intersect_point_1, intersect_point_2]) and
                                lib.is_points_on_face(target_fmw_face, [intersect_point_1, intersect_point_2])):
                                    intersected_points.append([intersect_point_1.X ,intersect_point_1.Y, intersect_point_1.Z])
                                    intersected_points.append([intersect_point_2.X ,intersect_point_2.Y, intersect_point_2.Z])
                if len(intersected_points) == 4:
                    wall_frw_obj.fmw_face.vertices = intersected_points
                    wall_frw_obj.fmw_face.point_min.X = round(intersected_points[0][0], 3)
                    wall_frw_obj.fmw_face.point_min.Y = round(intersected_points[0][1], 3)
                    wall_frw_obj.fmw_face.point_min.Z = round(intersected_points[0][2], 3)
                    wall_frw_obj.fmw_face.point_max.X = round(intersected_points[3][0], 3)
                    wall_frw_obj.fmw_face.point_max.Y = round(intersected_points[3][1], 3)
                    wall_frw_obj.fmw_face.point_max.Z = round(intersected_points[3][2], 3)
                    wall_frw_obj.fmw_face.ifc_space = self
                    wall_frw_obj.fmw_face.area = lib.polygon_area(intersected_points[0],intersected_points[3],abs(target_fmw_face.point_max.Z - target_fmw_face.point_min.Z))
                self.fmw_walls.append(wall_frw_obj)         
    def find_ifc_neighbour_spaces(self):
        # print(len(list(wall.ProvidesBoundaries)))
        # get all related spaces
        for fmw_wall in self.fmw_walls:
            for pb in fmw_wall.ifc_wall.ProvidesBoundaries: # rbe --> ifcRelSpaceBoundary
                if not (pb.RelatingSpace in self.neighbour_ifc_spaces) and (pb.RelatingSpace != self.ifc_space):
                    self.neighbour_ifc_spaces.append(pb.RelatingSpace)
    def find_fmw_adjacent_spaces(self):
        # floor_wall_dics = lib.find_floor_wall_faces(1)
        for target_fmw_wall in self.fmw_walls: # include an internal face
            # include an all faces of wall
            if target_fmw_wall.guid in config.floor_wall_dics.keys():
                individual_wall_faces_lst = config.floor_wall_dics[target_fmw_wall.guid] #{key: wall, values: faces} 
                p1 = target_fmw_wall.fmw_face.point_min
                p2 = target_fmw_wall.fmw_face.point_max
                main_line = LineString([(p1.X, p1.Y), (p2.X, p2.Y)])
                for f in individual_wall_faces_lst:
                    segment_points = [] # add the points of segmnets if two lines have common segment
                    if f.name != target_fmw_wall.fmw_face.name:
                        ## project the other faces line into target face line
                        p3 = lib.project_point_to_line(SHPoint(f.point_min.X, f.point_min.Y), main_line)
                        p4 = lib.project_point_to_line(SHPoint(f.point_max.X, f.point_max.Y), main_line)
                        p1_between = False
                        p2_between = False
                        p3_between = (min(p1.X, p2.X) <= round(p3[0], 3) <= max(p1.X, p2.X)) and (min(p1.Y, p2.Y) <= round(p3[1],3) <= max(p1.Y, p2.Y))
                        if p3_between:
                            segment_points.append(p3)
                        p4_between = (min(p1.X, p2.X) <= round(p4[0], 3) <= max(p1.X, p2.X)) and (min(p1.Y, p2.Y) <= round(p4[1],3) <= max(p1.Y, p2.Y))
                        if p4_between:
                            segment_points.append(p4)
                        if len(segment_points) < 2:
                            p1_between = (min(round(p3[0], 3), round(p4[0], 3)) <= p1.X <= max(round(p3[0], 3), round(p4[0], 3))) and (min(round(p3[1], 3), round(p4[1], 3)) <= p1.Y <= max(round(p3[1], 3), round(p4[1], 3)))
                            if p1_between:
                                segment_points.append([p1.X, p1.Y])
                            p2_between = (min(round(p3[0], 3), round(p4[0], 3)) <= p2.X <= max(round(p3[0], 3), round(p4[0], 3))) and (min(round(p3[1], 3), round(p4[1], 3)) <= p2.Y <= max(round(p3[1], 3), round(p4[1], 3)))
                            if p2_between:
                                segment_points.append([p2.X, p2.Y])
                        # if there is any intersection
                        # if p1_between or p2_between or p3_between or p4_between:
                        #     if f.ifc_space not in self.neighbour_fmw_spaces:
                        #         self.neighbour_fmw_spaces.append(f.ifc_space)
                        distance = 0
                        common_area = 0
                        if len(segment_points) == 2:
                            px = segment_points[0][0]
                            py = segment_points[0][1]
                            qx = segment_points[1][0]
                            qy = segment_points[1][1]
                            distance = round(math.sqrt(((px-qx)**2)+((py-qy)**2)),3)
                            common_area = round((distance * abs(p1.Z - p2.Z)), 3)
                            if f.ifc_space is not None and f.ifc_space.name in self.neighbour_fmw_spaces.keys():
                                self.neighbour_fmw_spaces[f.ifc_space.name] = self.neighbour_fmw_spaces[f.ifc_space.name] + common_area
                            elif f.ifc_space is not None:
                                self.neighbour_fmw_spaces[f.ifc_space.name] = common_area
class Wall(object):
    guid = ''
    tag = ''
    name = ''
    long_name = ''
    fmw_face = None
    area = 0
    def __init__(self, ifc_wall):
        self.ifc_wall = ifc_wall
        self.is_external = self.if_external()
        self.guid = self.ifc_wall.GlobalId
    def if_external(self):
        wall_property_set = None
        for definition in self.ifc_wall.IsDefinedBy:
            # To support IFC2X3, we need to filter our results.
            if definition.is_a('IfcRelDefinesByProperties'):
                property_set = definition.RelatingPropertyDefinition
                if(property_set.Name == 'Pset_WallCommon'): # Might return Pset_WallCommon
                    wall_property_set = definition.RelatingPropertyDefinition
                    break
            if wall_property_set != None:
                for property in wall_property_set.HasProperties:
                    if property.is_a('IfcPropertySingleValue') and property.Name == 'IsExternal':
                        return property.NominalValue.wrappedValue
class Face(object):
    wall = None
    name = None
    bbox = None
    area = 0
    point_min = None
    point_max = None
    normal_vector = None
    ifc_space = None
    ifc_face = None
    fmw_face = None
    vertices = None
class Point(object):
    def __init__(self, x, y, z):
        '''Defines x, y and z variables'''
        self.X = x
        self.Y = y
        self.Z = z
class Polygon(object):
    """ Polygon constructed from greater than two points.
    
        Only convex polygons are allowed! 
    
        Order of points is of course important!
    """
    ''' sample: points = [[0,0,0],[0,0.1,0],[0.1,0.1,-0.03],[0.1,0,-0.03]] '''
    def __init__(self,points):
        super(Polygon, self).__init__()
        self.pts = points
        #check if points are in one plane
        assert len(self.pts) >= 3, "You need at least 3 points to build a Polygon"
        if len(self.pts) > 3:
            x_0 = np.array(self.pts[0])
            for i in range(1,len(self.pts)-2):
                #the determinant of the vectors (volume) must always be 0
                x_i = np.array(self.pts[i])
                x_i1 = np.array(self.pts[i+1])
                x_i2 = np.array(self.pts[i+2])
                det = np.linalg.det([x_0-x_i, x_0-x_i1, x_0-x_i2])
                assert lib.cmp_floats( det, 0.0 ), "Points must be in a plane to create a Polygon"
                
    def on_surface(self, point):
        """Returns True if the point is on the polygon's surface and false otherwise."""
        n = len(self.pts)
        anglesum = 0
        p = np.array(point)

        for i in range(n):
            v1 = np.array(self.pts[i]) - p
            v2 = np.array(self.pts[(i+1)%n]) - p

            m1 = lib.magnitude(v1)
            m2 = lib.magnitude(v2)

            if lib.cmp_floats( m1*m2 , 0. ):
                return True #point is one of the nodes
            else:
                # angle(normal, vector)
                costheta = np.dot(v1,v2)/(m1*m2)
            anglesum = anglesum + np.arccos(costheta)
        return lib.cmp_floats( anglesum , 2*np.pi )


    def contains(self, point):
        return False


    def surface_identifier(self, surface_point, assert_on_surface = True):
        return "polygon"


    def surface_normal(self, ray, acute=False):
        vec1 = np.array(self.pts[0])-np.array(self.pts[1])
        vec2 = np.array(self.pts[0])-np.array(self.pts[2])
        normal = lib.norm(np.cross(vec1,vec2) )
        return normal
    def intersection(self, ray):
        """Returns a intersection point with a ray and the polygon."""
        n = self.surface_normal(ray)

        #Ray is parallel to the polygon
        if lib.cmp_floats( np.dot( np.array(ray.direction), n ), 0. ):
            return None
 
        t = 1/(np.dot(np.array(ray.direction),n)) * ( np.dot(n,np.array(self.pts[0])) - np.dot(n,np.array(ray.position)) )
        
        #Intersection point is behind the ray
        if t < 0.0:
            return None

        #Calculate intersection point
        point = np.array(ray.position) + t*np.array(ray.direction)
        
        #Check if intersection point is really in the polygon or only on the (infinite) plane
        if self.on_surface(point):
            return [list(point)]

        return None