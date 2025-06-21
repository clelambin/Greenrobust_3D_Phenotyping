"""Utilitiy functions used for the plant_volume modules"""

# Import libraries
import bpy               # Blender python
import bmesh             # Blender mesh module
import mathutils         # Blender object type
import numpy as np       # Array and matrix operations
from queue import Queue  # Multi-threading

# Complex type anotation
BMEntry = bmesh.types.BMVert | bmesh.types.BMEdge | bmesh.types.BMFace

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

def make_active(obj:bpy.types.Object) -> None:
    """Make input object as active and the only selected one"""
    bpy.ops.object.select_all(action='DESELECT')
    bpy.context.view_layer.objects.active = obj
    obj.select_set(True)

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
def get_connected_elements(input, cluster_type:str="faces", link_type:str="edges"):
    """Return a set containing the elements connected to the input element"""
    # Map cluster type to corresponding link object
    # eg : faces -> link_faces
    #      edges -> link_edges
    linked_input = f"link_{cluster_type}"
    return { c for l in getattr(input,link_type) for c in getattr(l, linked_input) if c != input }

# From https://blenderartists.org/t/get-amount-of-connected-geometry-within-a-mesh/1454143/2
def get_biggest_cluster(mesh:bmesh.types.BMesh, cluster_type:str="faces"):
    """Return biggest set of connected faces from input mesh"""
    # define link element used to connect cluster type entries
    link_entry = {"faces": "edges",
                  "edges": "verts"}
    link_type = link_entry[cluster_type]
    # Initialise mesh entry
    connected_groups = []
    work_list = [f for f in getattr(mesh, cluster_type)]
    while work_list:
        # Create worker for parrallel computation
        frontier = Queue()
        frontier.put( work_list[0] )
        this_group = [work_list[0]]
        work_list.pop(0)
        # Loop through all connected elements set to group them into one
        while not frontier.empty():
            connected_elements = get_connected_elements(frontier.get(), cluster_type, link_type)
            for next_element in connected_elements:
                if next_element not in this_group:
                    frontier.put(next_element)
                    this_group.append(next_element)
                    work_list.remove(next_element)
        # Append cluster to the list
        connected_groups.append(this_group)
    # Sort connected group by number of elements
    connected_groups = sorted(connected_groups, key=len)
    # Return biggest set of connected elements
    return connected_groups[-1]

def keep_biggest_cluster(mesh:bmesh.types.BMesh, cluster_type:str="faces") -> None:
    """Extract biggest mesh cluster from given object and delete vertices not part of the cluster"""
    # Extract biggest face cluster from given mesh
    cluster = get_biggest_cluster(mesh, cluster_type)
    # Extract set of all vertices included in mesh cluster
    # (Use double list comprehension: element.verts return a containing containing both vertices)
    vert_cluster = { vertex for element in cluster for vertex in element.verts }
    vert_all     = set(mesh.verts)
    # Use difference to get non included vertices
    vert_outside = vert_all.difference(vert_cluster)
    # Delete verts outside of face cluster
    delete_vertices(mesh, list(vert_outside))
