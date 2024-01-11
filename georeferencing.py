import numpy as N 

def augment(a):
    """Add a final column of ones to input data"""
    arr = N.ones((a.shape[0],a.shape[1]+1))
    arr[:,:-1] = a
    return arr

class Affine(object):
    def __init__(self, array=None):
        self.trans_matrix = array

    def transform(self, points):
        """Transform locally projected data using transformation matrix"""
        return N.dot(augment(N.array(points)), self.trans_matrix)

    @classmethod
    def from_tiepoints(cls, fromCoords, toCoords):
        "Produce affine transform by ingesting local and georeferenced coordinates for tie points"""
        fromCoords = augment(N.array(fromCoords))
        toCoords = N.array(toCoords)
        trans_matrix, residuals, rank, sv = N.linalg.lstsq(fromCoords, toCoords)

        affine =  cls(trans_matrix) # Setup affine transform from transformation matrix
        sol = N.dot(fromCoords,affine.trans_matrix) # Compute model solution
        print("Pixel errors:") 
        print (toCoords - sol)
        return affine
if __name__ == '__main__':
    input_file = r'C:\Users\Malavan-PC\OneDrive - York University\Data\Trapelo\IFC\space_xyz_with_att.txt'
    output_file = r'C:\Users\Malavan-PC\OneDrive - York University\Data\Trapelo\IFC\space_xyz_with_att_georef.txt'
    
    gps_points_local = [[219893.6501077091, 907158.5175406376],
                        [219915.51735378604, 907166.3054051218],
                        [219932.40876629358, 907191.9586738085],
                        [219913.0624046544, 907202.2576210345]]
    gps_points_geo = [[-71.2584331429999,42.414634253],
                        [-71.25809624,42.4146627650001],
                        [-71.25785114,42.415011038],
                        [-71.2581399329999, 42.41512593]]
    
    local_points_file = open(input_file, 'r')
    all_lines = local_points_file.readlines()
    local_points = []
    for line in all_lines:
        coords =[]
        coords.append(line.split(',')[0])
        coords.append(line.split(',')[1])
        local_points.append(coords)
    ### Transformation
    transform = Affine.from_tiepoints(gps_points_local,gps_points_geo)
    projected_data = transform.transform(local_points)
    i = 0
    with open(output_file, 'a') as out_file:
        for coor in projected_data:
            out_file.writelines("%s\n" % (str(coor[0]) + ',' + str(coor[1]) + ','+ all_lines[i].split(',')[2].rstrip() + ','+  all_lines[i].split(',')[3].rstrip()))
            i += 1
    out_file.close()