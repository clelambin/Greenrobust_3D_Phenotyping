"""Compute metrics relevant to the plant geometry
Measured metrics: - Volume
                  - Total surface area
                  - Dimension (X, Y and Z length of bounding box)
"""

# Import libraries
import bpy            # Blender python
import bmesh          # Blender mesh module
import numpy as np    # Array and matrix operations

# Import user modules
from cc_blender_plant_volume import blender_utility_functions as utility

def calc_area(mesh:bmesh.types.BMesh) -> float:
    """Return sum of all faces area"""
    return sum([face.calc_area() for face in mesh.faces])

def calc_dimension(mesh: bmesh.types.BMesh) -> np.ndarray:
    """Return X, Y and Z dimension of the bounding box of input mesh"""
    # Compute min and max for all 3 dimenions of the vertices coordinate
    vertex_coord = np.array([vertex.co for vertex in mesh.verts])
    return np.max(vertex_coord, axis=0) - np.min(vertex_coord, axis=0)

# Warning: if switch scaling to bmesh, would need to update
def calc_metrics(obj:bpy.types.Object) -> dict:
    """Return metrics (volume, area, dimension) from bmesh after application of transformation
    Return metrics are saved within a dictionary
    """
    # Mark object as active and switch to Edit mode
    bpy.ops.object.mode_set(mode='OBJECT', toggle=False)
    utility.select_all(select=False)
    bpy.context.view_layer.objects.active = obj
    bpy.ops.object.mode_set(mode='EDIT', toggle=False)
    # Extract bmesh from object
    orig_mesh = bmesh.from_edit_mesh(obj.data)
    temp_mesh = orig_mesh.copy()
    # Apply mesh transformation
    temp_mesh.transform(obj.matrix_world)
    # Calculate metrics
    volume = temp_mesh.calc_volume()
    surface = calc_area(temp_mesh)
    dimension = calc_dimension(temp_mesh)
    # Free up bmesh and switch back to Object mode
    orig_mesh.free()
    temp_mesh.free()
    bpy.ops.object.mode_set(mode='OBJECT', toggle=False)
    # Return all computed metrix
    return {"Volume"  : volume,
            "Surface" : surface,
            "Dim_X"   : dimension[0],
            "Dim_Y"   : dimension[1],
            "Dim_Z"   : dimension[2],}

if __name__ == "__main__":
    obj = bpy.context.active_object
    assert obj is not None, "No active object"
    plant_metrics = calc_metrics(obj)
    print(f"{plant_metrics = }")
