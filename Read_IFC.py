import sys
import ifcopenshell
def extract_quantities(property_definition):
    area = 0
    volume = 0
    length = 0
    print(property_definition.is_a())
    if 'IfcElementQuantity' == property_definition.is_a():
        for quantity in property_definition.Quantities:
            print(quantity)
            if 'IfcQuantityArea' == quantity.is_a():
                area =   quantity.AreaValue
            if 'IfcQuantityVolume' == quantity.is_a():
                volume = quantity.VolumeValue
            if 'IfcQuantityLength' == quantity.is_a():
                length = quantity.LengthValue
    return(area,volume,length)

if __name__ == "__main__":
    
    path = r'C:\Users\Malavan-PC\OneDrive - York University\Data\IFC\DuplexModel-IFC\Duplex_A_20110505.ifc'
    # path = r'C:\Users\Malavan-PC\OneDrive - York University\Data\Trapelo\IFC\Trapelo_Design.ifc'
    # path = r'C:\Users\Malavan-PC\OnseDrive - York University\Data\IFC\Samples\Duplex_A_20110505.ifc'
    ifc_file = ifcopenshell.open(path)
    # spaces = ifc_file.by_type('IfcSpace')
    space = ifc_file['0BTBFw6f90Nfh9rP1dlXrc']
    rbe = (space.BoundedBy)[10]
    element = rbe.RelatedBuildingElement
    # print(element.Name)
    print(element.is_a())
    # print(element.get_info())
    # for rbe in space.BoundedBy:
    #     element = rbe.RelatedBuildingElement
        # print(element.is_a())
        # print(element.Name)
        # print(element.get_info())
        # print(rbe.ConnectionGeometry)
    # floors = ifc_file.by_type("IfcBuildingStorey")
    # element = ifc_file['0QAEKuSavFZgLKbAPbK0ZW']
    # print(element)
    # space_boundaries = ifc_file.by_type('IfcRelSpaceBoundary')
    # for item in space_boundaries:
    #     print(item.name)
    # print(space.get_info())
    # products[0] ==  ifc_file[122](id) == ifc_file['2XQ$n5SLP5MBLyL442paFx'] (GlobalId)
    # obj_info = product.get_info()
    # print(obj_info.items())
    # for element in products:
    
    if element.IsDefinedBy:
        definitions = element.IsDefinedBy
        for definition in definitions:
            if 'IfcRelDefinesByProperties' == definition.is_a():
                property_definition = definition.RelatingPropertyDefinition
                # if 'IfcPropertySet' == property_definition.is_a():
                    # print(property_definition)
                geo_Props = extract_quantities(property_definition)
                elem_Props = element.get_info()
                # print('Id: {} -Name = {} - Area = {} - Volume = {} - Height = {}'.format(elem_Props['id'],elem_Props['Name'],geo_Props[0],geo_Props[1],geo_Props[2]))
            if 'IfcRelDefinesByType' == definition.is_a():
                type = definition.RelatingType
                # print(type)
                if type.HasPropertySets:
                    for property_definition in type.HasPropertySets:
                        extract_quantities(property_definition)