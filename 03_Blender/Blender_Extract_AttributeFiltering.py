# AttributeFiltering: Perform attribute based filtering on model (confidence, color, ...)

# Import libraries
from collections.abc import Callable
from queue import Queue  # Multi-threading
import bpy               # Blender python
import bmesh             # Blender mesh module
import mathutils         # Blender object type

# Complex type anotation
BMEntry = bmesh.types.BMVert | bmesh.types.BMEdge | bmesh.types.BMFace

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
    verts_thresh = [vertex for vertex in mesh.verts
                    if is_in_range(vertex[mesh_attr], min_, max_, transf)]
    return verts_thresh

def select_mesh_entry(mesh_list:list[BMEntry]) -> None:
    """Select vertices, edges or faces from mesh_list"""
    # Deselect all
    bpy.ops.mesh.select_all(action='DESELECT')
    # Select all mesh elements within the list
    for element in mesh_list:
        element.select_set(True)

def delete_vertices(mesh: bmesh.types.BMesh,
                    vertex_list:list[bmesh.types.BMVert]) -> None:
    """Delete vertices from vertex list"""
    for vertex in vertex_list:
        mesh.verts.remove(vertex)

# From https://blenderartists.org/t/get-amount-of-connected-geometry-within-a-mesh/1454143/2
def get_connected_faces(face):
    """Return a set containing the faces connected to the input face"""
    return { f for e in face.edges for f in e.link_faces if f != face }

# From https://blenderartists.org/t/get-amount-of-connected-geometry-within-a-mesh/1454143/2
def get_biggest_cluster(mesh:bmesh.types.BMesh):
    """Return biggest set of connected faces from input mesh"""
    connected_groups = []
    work_list = [f for f in mesh.faces]
    while work_list:
        # Create worker for parrallel computation
        frontier = Queue()
        frontier.put( work_list[0] )
        this_group = [work_list[0]]
        work_list.pop(0)
        # Loop through all connected faces set to group them into one
        while not frontier.empty():
            for next_face in get_connected_faces(frontier.get()):
                if next_face not in this_group:
                    frontier.put(next_face)
                    this_group.append(next_face)
                    work_list.remove(next_face)
        # Append face cluster to the list
        connected_groups.append(this_group)
    # Sort connected group by number of faces
    connected_groups = sorted(connected_groups, key=len)
    # Return biggest set of connected faces
    return connected_groups[-1]

def keep_biggest_cluster(mesh:bmesh.types.BMesh) -> None:
    """Extract biggest mesh cluster from given object and delete vertices not part of the cluster"""
    # Extract biggest face cluster from given mesh
    face_cluster = get_biggest_cluster(mesh)
    # Extract set of all vertices included in mesh cluster
    # (Use double list comprehension: face.verts return a containing containing both vertices)
    vert_cluster = { vertex for face in face_cluster for vertex in face.verts }
    vert_all     = set(mesh.verts)
    # Use difference to get non included vertices
    vert_outside = vert_all.difference(vert_cluster)
    # Delete verts outside of face cluster
    delete_vertices(mesh, list(vert_outside))

def extract_component(name:str="Pot",
                      color_selection:Callable=rgba_to_grey,
                      color_min:float=0,
                      color_max:float=float("Inf")) -> bpy.types.Object:
    """Duplicate active object and keep only selected component (selected by color range)"""
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
                                          min_=color_min,
                                          max_=color_max,
                                          layer="float_color",
                                          transf=color_selection)
    # Delete these vertices
    delete_vertices(obj_mesh, vertices_nonpot)
    # From remaining vertices, keep biggest face cluster
    keep_biggest_cluster(obj_mesh)
    # Switch back to object mode
    bpy.ops.object.mode_set(mode='OBJECT', toggle=False)
    # Return the extracted object
    return obj

def copy_pot() -> tuple[bpy.types.Object, bpy.types.Object]:
    """Extract pot and cup from active object (selected by color range)"""
    # Keep active object in memory (for second duplication)
    source_obj = bpy.context.active_object
    # Extract pot
    pot = extract_component(name="Pot", color_selection=rgba_to_grey, color_min=0.1)
    # Marke original active object as active and perform second copy
    bpy.ops.object.select_all(action='DESELECT')
    bpy.context.view_layer.objects.active = source_obj
    source_obj.select_set(True)
    cup = extract_component(name="Cup", color_selection=filter_pure_red, color_min=-0.1)
    # Return list of generated object
    return (pot, cup)

if __name__ == "__main__":
    # Test pot detection
    copy_pot()

#    # Test vertex selection tool
#    # Test on active object
#    obj = bpy.context.active_object
#    # Switch to Edit mode and load mesh
#    bpy.ops.object.mode_set(mode='EDIT', toggle=False)
#    obj_mesh = bmesh.from_edit_mesh(obj.data)
#    # Select biggest cluster of connected faces among selected mesh
#    face_group = get_biggest_cluster(obj_mesh)
#    select_mesh_entry(face_group)
#    print(f"{type(face_group)=}")
