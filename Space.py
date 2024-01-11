import timeit
import ifcopenshell
import IFC.Entities as model
import IFC.framework_library as lib
from shapely.geometry import Polygon
import IFC.config as config
if __name__ == '__main__':
    print('  -- Start the process...')
    start = timeit.default_timer() 
    path = r'C:\Users\Malavan-PC\OneDrive - York University\Data\IFC\DuplexModel-IFC\Duplex_A_20110505.ifc'
    # path = r'C:\Users\Malavan-PC\OneDrive - York University\Data\Trapelo\IFC\Trapelo_Design.ifc'
    ifc_file = ifcopenshell.open(path)
    all_ifc_spaces = ifc_file.by_type("IfcSpace")
    ifc_space = ifc_file['0BTBFw6f90Nfh9rP1dlXrc']
    floor_number = 2
    config.floor_wall_dics = lib.find_floor_wall_faces(floor_number,all_ifc_spaces)
    fmw_space = model.Space(ifc_space)
    # fmw_space.find_fmw_neighbour_spaces(floor_number)

    # for k, v in fmw_space.neighbour_fmw_spaces.items():
    #     print(k, v)

    for ifc_space in all_ifc_spaces: # all spaces (rooms) in the second floor
        if (ifc_space.Decomposes[0][4][2]).lower() == 'floor ' + str(floor_number):
            fmw_space = model.Space(ifc_space)
            fmw_space.find_fmw_adjacent_spaces()
            print('Framework - space: %s - neighbours: ' % ifc_space.Name)

            for k, v in fmw_space.neighbour_fmw_spaces.items():
                print(k, v)
    # print('IFC space neighbours: ' + ifc_space.Name)
    # for s in fmw_space.neighbour_ifc_spaces:
    #     print('  --%s' %s.Name)
    
    end = timeit.default_timer()
    print((' -- End of process. Duration:', round(end -start,4)))
    


