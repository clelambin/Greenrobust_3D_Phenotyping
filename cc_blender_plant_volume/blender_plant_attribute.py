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
import bpy            # Blender python
import bmesh          # Blender mesh module
import mathutils      # Blender object type

# Import user modules
from cc_blender_plant_volume import blender_point_clustering as cluster
from cc_blender_plant_volume import blender_extract_attributefiltering as attribute

# TODO: use decorator to extract mesh
def remove_ground(obj:bpy.types.Object, ground_height:float=0.002) -> None:
    """Remove all points from mesh below ground height"""
    # Extract mesh from object
    mesh = bmesh.new()
    mesh.from_mesh(obj.data)
    # Exract list of vertices with z < ground_height
    for vertex in mesh.verts:
        if vertex.co[2] < ground_height:
            mesh.verts.remove(vertex)
    # Apply modification to object
    mesh.to_mesh(obj.data)
    mesh.free()

def section_modifier(to_intersect:bpy.types.Object,
                     intersection_plane:bpy.types.Object,
                     thresh:float=0.0001) -> None:
    """Create intersection between object to intersect and plane
    Using boolean mesh operator
    """
    # Boolean modifier using fast intersection mode
    modifier = intersection_plane.modifiers.new(name="Section", type="BOOLEAN")
    modifier.solver = "FAST"
    modifier.operation = "INTERSECT"
    modifier.double_threshold = thresh
    # Set object to intersect
    modifier.object = to_intersect
    # Apply modifier
    bpy.ops.object.modifier_apply(modifier=modifier.name)

def create_plane(normal_axis:Literal["X", "Y", "Z"]="Z",
                 axis_position:float=0) -> bpy.types.Object:
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
    bpy.ops.mesh.primitive_plane_add(size=10, align='WORLD', location=axis_position*unit_vect[normal_axis])
    plane = bpy.context.active_object
    # Return error if plane could not be created
    assert plane is not None, "Error, section plane not created"
    # Rotate plane to be normal to specified axis
    plane.rotation_euler.rotate(rotation_mat[normal_axis])
    # Return plane object
    return plane

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
        section_modifier(obj, plane)

def plant_cleanup(obj:bpy.types.Object) -> None:
    """Extract green plant from 3D model, remove support and ouput model attribute
    Model attribute: plant volume, plant height, area projection
    """
    # Remove vertices below ground and extract green plants
    remove_ground(obj)
    plant_green = attribute.copy_green_plant(obj, color_thresh=0.02)
    cluster.dbscan_filter(plant_green, dbscan_eps=0.001)

def main() -> None:
    """Main function, cleanup selected object and output attributes"""
    plant = bpy.context.active_object
    # If no active object alert the user
    assert plant is not None, "No active object"

    # Cleanup plant model
    plant_cleanup(plant)

if __name__ == "__main__":
    # Test model preparation on active file
    main()
