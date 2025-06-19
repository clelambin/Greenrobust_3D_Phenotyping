# AttributeFiltering: Perform attribute based filtering on model (confidence, color, ...)

# Import libraries
from collections.abc import Callable
import bpy               # Blender python
import bmesh             # Blender mesh module
import mathutils         # Blender object type

# Import user modules
from cc_blender_plant_volume import blender_utility_functions as utility

def rgba_to_grey(RGBA:mathutils.Vector) -> float:
    """Convert color from RGBA (values from 0 to 1 for each channel) to greyscale"""
    # Convertion from RGB to greyscale using the Colorimetric conversion to grayscale
    # source: https://en.wikipedia.org/wiki/Grayscale
    return 0.2125*RGBA[0] + 0.7154*RGBA[1] + 0.0721*RGBA[2]

def filter_pure_red(RGBA:mathutils.Vector) -> float:
    """Convert color from RGBA (values from 0 to 1 for each channel) to filter out high red and low green and blue"""
    return -2*RGBA[0] + 4*RGBA[1] + 4*RGBA[2]

def is_in_range(value:float, min_:float, max_:float, transf:(None|Callable)=None) -> bool:
    """Return true if value is within [min_, max_]"""
    # (min_ and max_ are named with trailing _ to differenciate from the min and max function))
    # If a transformation function is specified, apply it to the value before checking if in range
    if transf is not None:
        value = transf(value)
    # Return true if value within [min_, max_]
    return value <= max_ and value >= min_

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
    verts_thresh = [vertex for vertex in mesh.verts
                    if is_in_range(vertex[mesh_attr], min_, max_, transf)]
    return verts_thresh

def edit_active_object(name:(str|None)=None) -> tuple[bpy.types.Object, bmesh.types.BMesh]:
    """Extract object and mesh from active object and swtich to edit mode"""
    obj = bpy.context.active_object
    # If no active object, alert the user
    assert obj is not None, "No active object"
    if name is not None:
        obj.name = name
    # Switch to Edit mode and load mesh
    bpy.ops.object.mode_set(mode='EDIT', toggle=False)
    obj_mesh = bmesh.from_edit_mesh(obj.data)
    # Return both object and mesh
    return (obj, obj_mesh)

def remesh_block_modifier(obj:bpy.types.Object, octree_depth:int=8) -> None:
    """Use block remesh modifier to simplify the geometry and remove unconnected"""
    # Add remesh modifier
    modifier = obj.modifiers.new(name="Remesh", type="REMESH")
    modifier.mode = "BLOCKS"
    modifier.scale = 0.9
    modifier.octree_depth = octree_depth
    # Make input object as active and switch to object mode to apply the modifier
    bpy.context.view_layer.objects.active = obj
    bpy.ops.object.mode_set(mode='OBJECT', toggle=False)
    # Apply modifier
    bpy.ops.object.modifier_apply(modifier="Remesh")

def extract_component(name:str="Pot",
                      color_selection:Callable=rgba_to_grey,
                      delete_min:float=0,
                      delete_max:float=float("Inf"),
                      filtering_mth:str|None=None,
                      remesh_octree:float|None=0.8) -> bpy.types.Object:
    """Duplicate active object and keep only selected component (selected by color range)"""
    # Duplicate active object
    bpy.ops.object.mode_set(mode='OBJECT', toggle=False)
    bpy.ops.object.duplicate()
    # Rename duplicated object and switch to edit mode
    (obj, obj_mesh) = edit_active_object(name=name)
    # Get list of non-pot vertices (with non-dark color)
    vertices_nonpot = vertex_by_attribute(obj_mesh, "Col",
                                          min_=delete_min,
                                          max_=delete_max,
                                          layer="float_color",
                                          transf=color_selection)
    # Delete these vertices
    utility.delete_vertices(obj_mesh, vertices_nonpot)
    # From remaining vertices, keep biggest face cluster
    if filtering_mth == "remesh":
        # Switch back to object mode to apply block remesh
        bpy.ops.object.mode_set(mode='OBJECT', toggle=False)
        remesh_block_modifier(obj, octree_depth=remesh_octree)
    elif filtering_mth == "cluster":
        # Stay in edit mode, then switch to object mode after face clustering
        utility.keep_biggest_cluster(obj_mesh)
        bpy.ops.object.mode_set(mode='OBJECT', toggle=False)

    # Return the extracted object
    return obj

def copy_pot() -> tuple[bpy.types.Object, bpy.types.Object]:
    """Extract pot and cup from active object (selected by color range)"""
    # Keep active object in memory (for second duplication)
    source_obj = bpy.context.active_object
    # Extract pot
    pot = extract_component(name="Pot",
                            color_selection=rgba_to_grey,
                            delete_min=0.1,
                            filtering_mth="cluster")
    # Marke original active object as active and perform second copy
    bpy.ops.object.select_all(action='DESELECT')
    bpy.context.view_layer.objects.active = source_obj
    source_obj.select_set(True)
    cup = extract_component(name="Cup",
                            color_selection=filter_pure_red,
                            delete_min=-0.1,
                            filtering_mth="remesh",
                            remesh_octree=8)
    # Return list of generated object
    return (pot, cup)

def copy_green_plant(source_obj:bpy.types.Object) -> bpy.types.Object:
    """Extract non-black part of the source model to get green parts"""
    # Mark plant as active and selected
    bpy.ops.object.select_all(action='DESELECT')
    bpy.context.view_layer.objects.active = source_obj
    source_obj.select_set(True)
    # Create copy of plant excluding pot
    plant_green = extract_component(name="Plant_green",
                                    color_selection=rgba_to_grey,
                                    delete_max=0.1)
    return plant_green

def delete_low_confidence() -> None:
    """Delete vertices with low confidence value from active object"""
    # Read active object and mesh
    _, mesh = edit_active_object()
    # Get list of vertices with low confidence value
    vertices_low_conf = vertex_by_attribute(mesh, max_=0)
    # Select low confidence vertices (Debug)
    utility.select_mesh_entry(vertices_low_conf)
    # Delete low confidence vertices
    utility.delete_vertices(mesh, vertices_low_conf)
    # Switch back to object mode
    bpy.ops.object.mode_set(mode='OBJECT', toggle=False)

if __name__ == "__main__":
#    # Test removing low confidence
#    delete_low_confidence()
#    # Test pot detection
#    copy_pot()
    obj = bpy.context.active_object
    copy_green_plant(obj)
