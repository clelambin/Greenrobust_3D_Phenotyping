"""Utilitiy functions used for the plant_volume modules"""

# Import libraries
import os                # File manager
from queue import Queue  # Multi-threading
import bpy               # Blender python
import bmesh             # Blender mesh module
from mathutils import Vector, Matrix    # Blender object type

# Import user modules
from cc_blender_plant_volume.blender_user_types import BMEntry, BMFaceList, Cartesian

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
    select_all(select=False)
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

def delete_selected(mesh: bmesh.types.BMesh,
                    selected_status:bool=True) -> bmesh.types.BMesh:
    """Delete selected or unselected vertices from input bmesh file, return updated bmesh"""
    for vertex in mesh.verts:
        if vertex.select == selected_status:
            mesh.verts.remove(vertex)
    return mesh

# From https://blenderartists.org/t/get-amount-of-connected-geometry-within-a-mesh/1454143/2
def get_connected_elements(input_element, cluster_type:str="faces", link_type:str="edges"):
    """Return a set containing the elements connected to the input element"""
    # Map cluster type to corresponding link object
    # eg : faces -> link_faces
    #      edges -> link_edges
    linked_input = f"link_{cluster_type}"
    return { c for l in getattr(input_element,link_type) for c in getattr(l, linked_input) if c != input_element }

# From https://blenderartists.org/t/get-amount-of-connected-geometry-within-a-mesh/1454143/2
def get_biggest_cluster(mesh:bmesh.types.BMesh, cluster_type:str="faces") -> BMFaceList:
    """Return biggest set of connected faces from input mesh"""
    # define link element used to connect cluster type entries
    link_entry = {"faces": "edges",
                  "edges": "verts"}
    link_type = link_entry[cluster_type]
    # Initialise mesh entry
    connected_groups = []
    work_list = list(getattr(mesh, cluster_type))
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
    # If no element found, alert the user and return None
    if len(connected_groups) == 0:
        print("Warning, no connected element found, return None")
        return None
    # If at least one group found return biggest set of connected elements
    return connected_groups[-1]


def keep_biggest_cluster(mesh:bmesh.types.BMesh, cluster_type:str="faces") -> int:
    """Extract biggest mesh cluster from given object and delete vertices not part of the cluster"""
    # Extract biggest face cluster from given mesh
    cluster = get_biggest_cluster(mesh, cluster_type)
    # If cluster return None, indicate that no element found, skip the element deletion
    if cluster is None:
        return 0
    # Extract set of all vertices included in mesh cluster
    # (Use double list comprehension: element.verts return a containing containing both vertices)
    vert_cluster = { vertex for element in cluster for vertex in element.verts }
    vert_all     = set(mesh.verts)
    # Use difference to get non included vertices
    vert_outside = vert_all.difference(vert_cluster)
    # Delete verts outside of face cluster
    delete_vertices(mesh, list(vert_outside))
    # Return the amount of deleted vertices
    return len(vert_outside)

def select_all(select:bool= True, target_type:str="MESH") -> None:
    """Select or deselect all object of the given type"""
    # Loop through existing object, if object of the given type, edit selection
    for obj in bpy.data.objects:
        if obj.type == target_type:
            obj.select_set(select)

def get_view3d() -> bpy.types.Area | None:
    """Return first 3d view in Blender environment"""
    for area in bpy.context.screen.areas:
        if area.type == "VIEW_3D":
            return area
    # If no area of the VIEW_3D type found, return None
    return None

def get_active_obj() -> bpy.types.Object:
    """Return active object, raise an error of no object is selected"""
    obj = bpy.context.active_object
    assert obj is not None, "No active object"
    return obj

def hide_object(obj:bpy.types.Object, hide_bool:bool=True):
    """Set viewport and render visibility (if boolean set to True, hide, otherwide, show)"""
    obj.hide_set(hide_bool)
    obj.hide_render = hide_bool
    # hide_viewport hide in all viewport, not only in the current one
#    obj.hide_viewport = False

def from_z_to_axis(obj:bpy.types.Object, axis:Cartesian) -> None:
    """Rotate object to allign to given axis (starting from Z axis)"""
    # Define rotation matrix for each axis
    rotation_mat = {
        "X": Matrix([[0,0,-1],[0,1,0],[1,0,0]]),
        "Y": Matrix([[1,0,0],[0,0,-1],[0,1,0]]),
        "Z": Matrix([[1,0,0],[0,1,0],[0,0,1]]),
    }
    obj.rotation_euler.rotate(rotation_mat[axis])

def create_curve(point_list:list[Vector], name:str="Curve") -> bpy.types.Curve:
    """Create a blender curve object with points from point list"""
    # Initialise the curve object
    curve_data = bpy.data.curves.new(name, "CURVE")
    curve_data.dimensions = "3D"
    # Add polyline to curve with a point at center of each node of straight branch
    curve_poly = curve_data.splines.new("POLY")
    curve_poly.points.add(len(point_list) - 1)
    for (index, coord) in enumerate(point_list):
        # Polyline coodinate is composed of 4 components: x, y, z (part of coord) and the weight
        curve_poly.points[index].co = coord.to_4d()
    return curve_data

def create_plane(normal_axis:Cartesian="Z",
                 axis_position:float=0,
                 name_digits:int=3,
                 hide:bool=False) -> bpy.types.Object:
    """Create a plane normal to given axis and passing by normal axis as input position"""
    # Define unit vetor used to define position for each axis
    unit_vect = {
        "X": Vector([1,0,0]),
        "Y": Vector([0,1,0]),
        "Z": Vector([0,0,1]),
    }
    # Create plane at given location
    # (by default the plane is alligned with z)
    location = axis_position*unit_vect[normal_axis]
    bpy.ops.mesh.primitive_plane_add(size=10, align='WORLD', location=location)
    plane = get_active_obj()
    # .{name_digit}f add digit to the number to improve object order
    plane.name = f"Plane_{normal_axis}-{axis_position:.{name_digits}f}"
    # Return error if plane could not be created
    assert plane is not None, "Error, section plane not created"
    # Rotate plane to be normal to specified axis
    from_z_to_axis(plane, normal_axis)
    # If requested, hide the created object
    if hide:
        hide_object(plane)
    # Return plane object
    return plane

def curve_to_mesh(curve:bpy.types.Curve, del_curve:bool=True) -> bpy.types.Object:
    """Convert input curve into an mesh object and delete input curve (if requested)"""
    # Create object and link to scene
    curve_obj  = bpy.data.objects.new(curve.name, curve)
    curve_mesh = bpy.data.meshes.new_from_object(curve_obj)
    mesh_obj   = bpy.data.objects.new(curve.name, curve_mesh)
    bpy.context.collection.objects.link(mesh_obj)

    # Cleanup intermediate data
    if del_curve:
        bpy.data.objects.remove(curve_obj)
        bpy.data.curves.remove(curve)

    # Output converted object
    return mesh_obj

def extract_face_from_verts(verts:list[bmesh.types.BMVert]) -> bmesh.types.BMFace | None:
    """Extract face shared by input vertex"""
    vertices_set = set(verts)
    # Loop through faces connected to first vertex
    for face in verts[0].link_faces:
        # If set of vertices in current face equal to input vertices, return face
        if set(face.verts) == vertices_set:
            return face
    # No face found with same set of vertices, return None
    return None

def offset_faces(mesh:bmesh.types.BMesh, dist:float=0.001) -> None:
    """Offset faces of mesh by input distance"""
    for face in mesh.faces:
        # Renormalise the face normal
        offset_dist = face.normal.normalized()*dist
        bmesh.ops.translate(mesh, vec=offset_dist, verts = list(face.verts))

def cleanup_env(obj_to_remove:list[str] = ["Cube",], type_to_remove:list[str] = ["MESH",]) -> None:
    """Remove non-relevant object from scene"""
    # Loop through object and delete object in the list
    for obj in bpy.data.objects:
        if obj.type in type_to_remove or obj.name in obj_to_remove:
            bpy.data.objects.remove(obj)

def import_file(filepath:str, file_ext:str="obj") -> bpy.types.Object:
    """Import obj or ply file and return as an object"""
    # Check if file exist and is obj
    if not os.path.isfile(filepath):
        raise OSError(f"File {filepath} not found")
    if not filepath.endswith(file_ext):
        raise OSError(f"File {filepath} is not an obj")

    # Import file (based on the file_ext)
    import_command = {"ply":"ply_import", "obj":"obj_import"}
    getattr(bpy.ops.wm, import_command[file_ext])(filepath=filepath,
                                                  forward_axis='X',
                                                  up_axis='Z')
    # Apply rotation to consider for the imported axis transformation
    bpy.ops.object.transform_apply(location=False, rotation=True, scale=False)

    # Imported object is marked aas active, return it
    return get_active_obj()

def save_blend(obj_path:str, output_path:str) -> None:
    """Save the prepared model as blend file in output file"""
    # Extract plant name
    plant_name = obj_path.split(os.sep)[-1]
    plant_name = plant_name.split(".")[0]
    plant_file = f"{plant_name}.blend"
    bpy.ops.wm.save_as_mainfile(filepath=os.path.join(output_path, plant_file))

def scene_statistics() -> dict[str, str]:
    """Return the scene statistics of the active object as a dictionary"""
    # Scene statistics outputed as a full string like such:
    # 'Scene Collection | Metashape_BR07IC_20251110_1658 | Verts:4 | Faces:4 | Tris:4 | Objects:1/1 | Duration: 00:10+10 (Frame 1/250) | Memory: 74.6 MiB | 4.4.3'
    stats_string = bpy.context.scene.statistics(bpy.context.view_layer)
    stats_dict   = {}
    # Split by " | "
    for stat in stats_string.split(" | "):
        # If element is a property, it contain a : which split the name and the value
        property = stat.split(":")
        if len(property) == 2:
            stats_dict[property[0]] = property[1]
    # Return the dictionary containing the scene properties
    return stats_dict
