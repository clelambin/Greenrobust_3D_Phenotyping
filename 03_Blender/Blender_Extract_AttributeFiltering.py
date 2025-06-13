# AttributeFiltering: Perform attribute based filtering on model (confidence, color, ...)

# Import libraries
from collections.abc import Callable
import bpy            # Blender python
import bmesh          # Blender mesh module
import mathutils      # Blender object type

def rgba_to_grey(RGBA:mathutils.Vector) -> float:
    """Convert color from RGBA (values from 0 to 1 for each channel) to greyscale"""
    return 0.2125*RGBA[0] + 0.7154*RGBA[1] + 0.0721*RGBA[2]

def is_in_range(value:float, min_:float, max_:float, transf:(None|Callable)=None) -> bool:
    """Return true if value is within [min_, max_]"""
    # (min_ and max_ are named with trailing _ to differenciate from the min and max function))
    # If a transformation function is specified, apply it to the value before checking if in range
    if transf is not None:
        value = transf(value)
    # Return true if value within [min_, max_]
    return value < max_ and value > min_

def vertex_by_attribute(mesh:bmesh.types.BMesh,
                        attribute:str="confidence",
                        min_:float=0.0,
                        max_:float=float("Inf"),
                        layer:str="float",
                        transf:None|Callable=None) -> list[bmesh.types.BMVert]:
    """Return vertex list which attribute exceed the threshold for the given attribute"""
    # Load attribute to use for filtering
    mesh_attr = getattr(mesh.verts.layers, layer)[attribute]
    # Return vertex whose attribute is bounded between min and max
    verts_thresh = [vertex for vertex in mesh.verts if is_in_range(vertex[mesh_attr], min_, max_, transf)]
    return verts_thresh

def select_vertices(mesh: bmesh.types.BMesh,
                    vertex_list:list[bmesh.types.BMVert]) -> None:
    """Select vertices from vertex_list"""
    # Deselect all
    bpy.ops.mesh.select_all(action='DESELECT')
    # Select all vertices within the list
    for vertex in vertex_list:
        vertex.select_set(True)

def delete_vertices(mesh: bmesh.types.BMesh,
                    vertex_list:list[bmesh.types.BMVert]) -> None:
    """Delete vertices from vertex list"""
    for vertex in vertex_list:
        mesh.verts.remove(vertex)

def copy_pot(name:str="Pot") -> None:
    """Duplicate active object and keep only pot (selected by color range)"""
    # Duplicate active object
    bpy.ops.object.mode_set(mode='OBJECT', toggle=False)
    bpy.ops.object.duplicate()
    # Rename duplicated object
    obj = bpy.context.active_object
    obj.name = name
    # Switch to Edit mode and load mesh
    bpy.ops.object.mode_set(mode='EDIT', toggle=False)
    obj_mesh = bmesh.from_edit_mesh(obj.data)
    # Get list of non-pot vertices (with non-dark color)
    vertices_nonpot = vertex_by_attribute(obj_mesh, "Col",
                                          min_=0.1,
                                          layer="float_color",
                                          transf=rgba_to_grey)
    # Delete these vertices
    delete_vertices(obj_mesh, vertices_nonpot)
    # Switch back to object mode
    bpy.ops.object.mode_set(mode='OBJECT', toggle=False)

if __name__ == "__main__":
    # Test pot detection
    copy_pot()

#    # Test vertex selection tool
#    # Test on active object
#    obj = bpy.context.active_object
#    # Switch to Edit mode and load mesh
#    bpy.ops.object.mode_set(mode='EDIT', toggle=False)
#    obj_mesh = bmesh.from_edit_mesh(obj.data)
#    # Select vertex with a too low confidence
#    vertices_low_conf = vertex_by_attribute(obj_mesh, "confidence", max_=10.0)
#    select_vertices(obj_mesh, vertices_low_conf)
#    # Select dark vertices (represent the pot)
#    vertices_pot = vertex_by_attribute(obj_mesh, "Col",
#                                       max_=0.1,
#                                       layer="float_color",
#                                       transf=rgba_to_grey)
#    select_vertices(obj_mesh, vertices_pot)
