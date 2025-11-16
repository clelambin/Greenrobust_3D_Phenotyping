"""Main script for the model processing from already alligned and scaled model

Note:
- require digitisation with marker for allignment
- replace model_prep fromblender_plant_modelprep.py

Workflow:
- Import obj
- Use color attribute to extract green plant
- Remove support (wood stick or torus) from model
- Extract model attribute (volume, ...)
"""

# Import libraries
from typing import Literal
from math import pi
import bpy            # Blender python
import bmesh          # Blender mesh module
import mathutils      # Blender object type

# Import user modules
from cc_blender_plant_volume import blender_point_clustering as cluster
from cc_blender_plant_volume import blender_extract_attributefiltering as attribute
from cc_blender_plant_volume import blender_utility_functions as utility

def remove_ground(obj:bpy.types.Object, ground_height:float=0.002) -> None:
    """Remove all points from mesh below ground height"""
    # Add cube object, shift it toward the ground
    cube_location = mathutils.Vector([0, 0, ground_height-1])
    bpy.ops.mesh.primitive_cube_add(size=2, location=cube_location)
    ground = utility.get_active_obj()
    ground.name = "Ground"
    # Apply difference boolean to remove anything inside the ground object
    boolean_modifier(source_obj=obj, target_obj=ground, operation="DIFFERENCE", thresh=0)
    # Hide Ground
    utility.hide_object(ground)

def boolean_modifier(source_obj:bpy.types.Object,
                     target_obj:bpy.types.Object,
                     thresh:float=0.00001,
                     operation:Literal["DIFFERENCE", "INTERSECT"]="INTERSECT",
                     apply:bool=True)-> None:
    """Create intersection between object to intersect and plane
    Using boolean mesh operator
    """
    # Mark source object as active
    utility.make_active(source_obj)
    # Boolean modifier using fast intersection mode
    modifier = source_obj.modifiers.new(name="Section", type="BOOLEAN")
    modifier.solver = "FAST"
    modifier.operation = operation
    modifier.double_threshold = thresh
    # Set object to intersect
    modifier.object = target_obj
    # Apply modifier
    if apply:
        bpy.ops.object.modifier_apply(modifier=modifier.name)

def create_plane(normal_axis:Literal["X", "Y", "Z"]="Z",
                 axis_position:float=0,
                 name_digits:int=3) -> bpy.types.Object:
    """Create a plane normal to given axis and passing by normal axis as input position"""
    # Define unit vetor used to define position for each axis
    unit_vect = {
        "X": mathutils.Vector([1,0,0]),
        "Y": mathutils.Vector([0,1,0]),
        "Z": mathutils.Vector([0,0,1]),
    }
    # Define rotation matrix for each axis
    rotation_mat = {
        "X": mathutils.Matrix([[0,0,1],[0,1,0],[-1,0,0]]),
        "Y": mathutils.Matrix([[1,0,0],[0,0,-1],[0,1,0]]),
        "Z": mathutils.Matrix([[1,0,0],[0,1,0],[0,0,1]]),
    }
    # Create plane at given location
    # (by default the plane is alligned with z)
    location = axis_position*unit_vect[normal_axis]
    bpy.ops.mesh.primitive_plane_add(size=10, align='WORLD', location=location)
    plane = utility.get_active_obj()
    # .{name_digit}f add digit to the number to improve object order
    plane.name = f"Plane_{normal_axis}_{axis_position:.{name_digits}f}"
    # Return error if plane could not be created
    assert plane is not None, "Error, section plane not created"
    # Rotate plane to be normal to specified axis
    plane.rotation_euler.rotate(rotation_mat[normal_axis])
    # Return plane object
    return plane

def calc_rondness(face:bmesh.types.BMFace) -> float:
    """Use the Polsby–Popper to measure the rondness (or compactness) of a plannar face"""
    area = face.calc_area()
    perimeter = face.calc_perimeter()
    return (4 * pi * area) / (perimeter**2)

def section_faces(obj:bpy.types.Object) -> list[dict]|None:
    """Remove vertices not connected to any faces from current object mesh,
    Return a list containing the roundness and center of each face in the section
    """
    # initialise face description list
    face_descr = []

    # Extract mesh from object
    mesh = bmesh.new()
    mesh.from_mesh(obj.data)

    # If object do not contain any faces, delete the object
    if len(mesh.faces) == 0:
        bpy.data.objects.remove(obj, do_unlink=True)
        return None

    # Save all vertices connected to faces as a set
    connected_verts = set()
    for face in mesh.faces:
        connected_verts.update(face.verts)
        # Save roundness and center of current face
        face_descr.append({
            "roundness": calc_rondness(face),
            "center": face.calc_center_median(),
            "area": face.calc_area()
        })

    # Get vertices not parts of connected vertices and delete them
    all_verts = set(mesh.verts)
    disconnected_verts = list(all_verts.difference(connected_verts))
    utility.delete_vertices(mesh, disconnected_verts)

    # Apply modification to object
    mesh.to_mesh(obj.data)
    mesh.free()

    # Return the face description list, sorted by roundness
    return sorted(face_descr, key=lambda face: face["roundness"], reverse=True)

# TODO: check that obj is updated before getting dimension
def multi_crosssection(obj:bpy.types.Object, z_delta:float=0.1) -> None:
    """Create crosssection of input object spaced by z_delta"""
    # Get Z dimention to define number of cross-sections
    obj_dim    = obj.dimensions
    nb_section = int(obj_dim[2] / z_delta) + 1
    for section_index in range(1, nb_section):
        section_z = section_index * z_delta
        # Create plane at given section_z height
        plane = create_plane("Z", section_z)
        # Apply section modifier to plane
        boolean_modifier(source_obj=plane, target_obj=obj)
        # Cleanup cross-section and save list of all faces within section
        section_faces(plane)

def plant_cleanup(obj:bpy.types.Object) -> None:
    """Extract green plant from 3D model, remove support and ouput model attribute
    Model attribute: plant volume, plant height, area projection
    """
    # Remove vertices below ground and extract green plants
    remove_ground(obj)
#    plant_green = attribute.copy_green_plant(obj, color_thresh=0.02)
#    cluster.dbscan_filter(plant_green, dbscan_eps=0.001)

def main() -> None:
    """Main function, cleanup selected object and output attributes"""
    plant = bpy.context.active_object
    # If no active object alert the user
    assert plant is not None, "No active object"

    # Cleanup plant model
    plant_cleanup(plant)

    # Plant section
    multi_crosssection(plant, z_delta=0.02)

if __name__ == "__main__":
    # Test model preparation on active file
    main()
