import numpy as np
import ifcopenshell
import IFC.framework_library as lib
from shapely.geometry import Point
from shapely.geometry import LineString
import IFC.framework_library as lib
import ifcopenshell
import ifcopenshell.geom

def cmp_floats(a,b, atol=1e-12):
    return abs(a-b) < atol


def magnitude(vector):
   return np.sqrt(np.dot(np.array(vector), np.array(vector)))


def norm(vector):
   return np.array(vector)/magnitude(np.array(vector))

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
    print(P) #0.2 1.
class Ray(object):
    """A ray in the global cartesian frame."""
    def __init__(self, position, direction):
        self.position = np.array(position)
        self.direction = norm(direction) 


class Polygon(object):
    """ Polygon constructed from greater than two points.
    
        Only convex polygons are allowed! 
    
        Order of points is of course important!
    """
    
    def __init__(self, points):
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
                assert cmp_floats( det, 0.0 ), "Points must be in a plane to create a Polygon"
                

    def on_surface(self, point):
        """Returns True if the point is on the polygon's surface and false otherwise."""
        n = len(self.pts)
        anglesum = 0
        p = np.array(point)

        for i in range(n):
            v1 = np.array(self.pts[i]) - p
            v2 = np.array(self.pts[(i+1)%n]) - p

            m1 = magnitude(v1)
            m2 = magnitude(v2)

            if cmp_floats( m1*m2 , 0. ):
                return True #point is one of the nodes
            else:
                # angle(normal, vector)
                costheta = np.dot(v1,v2)/(m1*m2)
            anglesum = anglesum + np.arccos(costheta)
        return cmp_floats( anglesum , 2*np.pi )


    def contains(self, point):
        return False


    def surface_identifier(self, surface_point, assert_on_surface = True):
        return "polygon"


    def surface_normal(self, ray, acute=False):
        vec1 = np.array(self.pts[0])-np.array(self.pts[1])
        vec2 = np.array(self.pts[0])-np.array(self.pts[2])
        normal = norm( np.cross(vec1,vec2) )
        return normal


    def intersection(self, ray):
        """Returns a intersection point with a ray and the polygon."""
        n = self.surface_normal(ray)

        #Ray is parallel to the polygon
        if cmp_floats( np.dot( np.array(ray.direction), n ), 0. ):
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


if __name__ == '__main__':
    # path = r'C:\Users\Malavan-PC\OneDrive - York University\Data\IFC\DuplexModel-IFC\Duplex_A_20110505.ifc'
    path = r'C:\Users\Malavan-PC\OneDrive - York University\Data\Trapelo\IFC\Trapelo_Design.ifc'
    space_3D_file = open(r'C:\Users\Malavan-PC\OneDrive - York University\Data\Trapelo\IFC\space_xyz_with_att.txt',"w")
    ifc_file = ifcopenshell.open(path)
    all_ifc_spaces = ifc_file.by_type("IfcSpace")

    settings = ifcopenshell.geom.settings()
    settings.set(settings.USE_PYTHON_OPENCASCADE, True)

    for ifc_space in all_ifc_spaces: 
        space_bbx_cntr = lib.space_bbx_center(ifc_space,settings).XYZ()
        space_3D_file.write(str(space_bbx_cntr.X())+ ','+str(space_bbx_cntr.Y())+','+str(space_bbx_cntr.Z())+','+ ifc_space.Name + "\n")
    space_3D_file.close() 
    print(all_ifc_spaces)
   
