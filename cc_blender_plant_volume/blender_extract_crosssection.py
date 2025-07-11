"""Extract cross section of object at a given location, used for model scaling"""

# To do:
# - Find way to rename input socket (instead of using Socket_2)

import bpy
import bmesh
import numpy as np

# Import user modules
from cc_blender_plant_volume import blender_utility_functions as utility

def geonode_init(name:str = "Geometry node") -> bpy.types.GeometryNodeTree:
    """Initialise geo-node object"""
    # Create Geometry node
    geonode  = bpy.data.node_groups.new(type = 'GeometryNodeTree', name = name)
    geonode.is_modifier = True
    # Create input and output socket (the input and output of each node)
    geonode.interface.new_socket(name = "Geometry", in_out='INPUT', socket_type = 'NodeSocketGeometry')
    geonode.interface.new_socket(name = "Geometry", in_out='OUTPUT', socket_type = 'NodeSocketGeometry')
    # Create group input and output nodes
    group_input  = geonode.nodes.new("NodeGroupInput")
    group_output = geonode.nodes.new("NodeGroupOutput")
    group_input.name  = "Group Input"
    group_output.name = "Group Output"
    # Return initialised geometry node
    return geonode

def cross_section_node(name:str="Geometry node",
                       boolean_intersect:bool=True) -> bpy.types.GeometryNodeTree:
    """Create geometry node on plane object and use it to generate cross-section modifier"""
    # Initialise cross-section geometry node
    geonode = geonode_init(name=name)
    group_input  = geonode.nodes["Group Input"]
    group_output = geonode.nodes["Group Output"]
    # Create additional input socket to reference external object
    geonode.interface.new_socket(name = "Object", in_out='INPUT', socket_type = 'NodeSocketObject')
    # Add reference to external object
    object_info = geonode.nodes.new("GeometryNodeObjectInfo")
    object_info.name = "Object Info"
    object_info.transform_space = 'RELATIVE'
    # Add boolean geometry node (to create the cross-section)
    mesh_boolean = geonode.nodes.new("GeometryNodeMeshBoolean")
    mesh_boolean.name = "Mesh Boolean" 
    mesh_boolean.operation = 'INTERSECT' if boolean_intersect else 'DIFFERENCE'
    mesh_boolean.solver = 'EXACT'
    # Create links between nodes
    geonode.links.new(group_input.outputs["Object"], object_info.inputs["Object"])
    if boolean_intersect:
        geonode.links.new(group_input.outputs["Geometry"], mesh_boolean.inputs["Mesh"])
        geonode.links.new(object_info.outputs["Geometry"], mesh_boolean.inputs["Mesh"])
        geonode.links.new(mesh_boolean.outputs["Mesh"], group_output.inputs["Geometry"])
    else:
        # If use difference boolean, add separate geometry node (to isolate the cross-section)
        # and link it to the other nodes
        geonode.links.new(group_input.outputs["Geometry"], mesh_boolean.inputs["Mesh 1"])
        geonode.links.new(object_info.outputs["Geometry"], mesh_boolean.inputs["Mesh 2"])
        separate_geometry = geonode.nodes.new("GeometryNodeSeparateGeometry")
        separate_geometry.name = "Separate Geometry"
        separate_geometry.domain = 'POINT'
        geonode.links.new(mesh_boolean.outputs["Mesh"], separate_geometry.inputs["Geometry"])
        geonode.links.new(mesh_boolean.outputs["Intersecting Edges"], separate_geometry.inputs["Selection"])
        geonode.links.new(separate_geometry.outputs["Selection"], group_output.inputs["Geometry"])
    # Return geometry node
    return geonode

def cross_section_modifier(plane:bpy.types.Object, intersect:bpy.types.Object, name:str="Geometry node") -> None:
    """Apply cross-section modifier on input object to intersect"""
    # Create modifier
    modifier = plane.modifiers.new(name=name, type='NODES')
    # Check if geometry node exist, if not, create it
    if bpy.data.node_groups.find(name) >= 0:
        geonode = bpy.data.node_groups[name]
    else:
        geonode = cross_section_node(name=name)
    # Assign geometry node to modifier
    modifier.node_group = geonode
    # Assign external geometry socket
    modifier["Socket_2"] = intersect
    # Apply modifier
    bpy.ops.object.modifier_apply(modifier=modifier.name)

def keep_intersection(obj: bpy.types.Object) -> None:
    """Keep only intersection plane after cross-section modifier
    After cross-section modifier, only the vertices in the intersection are selected
    Use bmesh to delete non-selected vertices
    """
    # Extract bmesh from obj
    obj_mesh = bmesh.new()
    obj_mesh.from_mesh(obj.data)
    # Delete unselected vertices
    obj_mesh = utility.delete_selected(obj_mesh, selected_status=False)
    # Update object mesh
    obj_mesh.to_mesh(obj.data)
    obj_mesh.free()

def main(intersect:bpy.types.Object) -> bpy.types.Object:
    """Main script, use geometry node to create cross section of intersect object"""
    # Create plane at given location
    bpy.ops.mesh.primitive_plane_add(size=10, align='WORLD', location=(0, 0, 0))
    plane = bpy.context.active_object
#    # Create geometry node on plane
#    cross_section_node(name="Cross-section")
    # Apply modifier on intersect object
    cross_section_modifier(plane, intersect, name="Cross-section")
    keep_intersection(plane)
    # Return generated plane
    return plane

# From https://stackoverflow.com/questions/13542855/algorithm-to-find-the-minimum-area-rectangle-for-given-points-in-order-to-comput
def minimum_bounding_rectangle(points):
    """
    Find the smallest bounding rectangle for a set of points.
    Returns a set of points representing the corners of the bounding box.

    :param points: an nx2 matrix of coordinates
    :rval: an nx2 matrix of coordinates
    """
    pi2 = np.pi/2.
    # calculate edge angles
    edges = np.zeros((len(points)-1, 2))
    edges = points[1:] - points[:-1]
    angles = np.zeros((len(edges)))
    angles = np.arctan2(edges[:, 1], edges[:, 0])
    angles = np.abs(np.mod(angles, pi2))
    angles = np.unique(angles)
    # find rotation matrices
    rotations = np.vstack([np.cos(angles),
                           np.cos(angles-pi2),
                           np.cos(angles+pi2),
                           np.cos(angles)]).T
    rotations = rotations.reshape((-1, 2, 2))

    # apply rotations to the hull
    rot_points = np.dot(rotations, points.T)
    # find the bounding points
    min_x = np.nanmin(rot_points[:, 0], axis=1)
    max_x = np.nanmax(rot_points[:, 0], axis=1)
    min_y = np.nanmin(rot_points[:, 1], axis=1)
    max_y = np.nanmax(rot_points[:, 1], axis=1)

    # find the box with the best area
    areas = (max_x - min_x) * (max_y - min_y)
    best_idx = np.argmin(areas)
    # return the best box
    x1 = max_x[best_idx]
    x2 = min_x[best_idx]
    y1 = max_y[best_idx]
    y2 = min_y[best_idx]
    r = rotations[best_idx]
    rval = np.zeros((4, 2))
    rval[0] = np.dot([x1, y2], r)
    rval[1] = np.dot([x2, y2], r)
    rval[2] = np.dot([x2, y1], r)
    rval[3] = np.dot([x1, y1], r)
    # Return bounding box coordinates
    return rval

def draw_polygon_2d(vertex_coord, z:float=0.0) -> None:
    """Draw polygon from vertices list on active edit mesh"""
    # Get active object mesh
    obj = bpy.context.active_object
    obj_mesh = bmesh.from_edit_mesh(obj.data)
    # Loop through vertex coordinates and add vertices
    vertex_list = []
    for coord in vertex_coord:
        vertex = obj_mesh.verts.new((coord[0], coord[1], z))
        vertex_list.append(vertex)
    # Loop through vertex list and add edges
    vertex_prev = vertex_list[0]
    vertex_init = vertex_prev
    for vertex in vertex_list[1:]:
        obj_mesh.edges.new([vertex_prev, vertex])
        vertex_prev = vertex
    obj_mesh.edges.new([vertex_prev, vertex_init])
    # Update mesh and free bmesh memory
    bmesh.update_edit_mesh(obj.data)
    obj_mesh.free()

def remove_long_edges(mesh:bmesh.types.BMesh, max_length: float=0.1) -> None:
    """Remove edges exceeding max_length"""
    for edge in mesh.edges:
        v1, v2 = [vertex.co for vertex in edge.verts]
        dist = np.linalg.norm(np.array(v2-v1))
        if dist > max_length:
            mesh.edges.remove(edge)

def get_closest_to_center(vertex_list:bmesh.types.BMVertSeq) -> bmesh.types.BMVert:
    """Return vertex closest to center"""
    smallest_dist  = float("Inf")
    closest_vertex = None
    for vertex in vertex_list:
        dist_to_center = np.linalg.norm(np.array(vertex.co))
        if dist_to_center < smallest_dist:
            smallest_dist = dist_to_center
            closest_vertex = vertex
    # If no closest vertex found, alert the user
    assert closest_vertex is not None, "No vertex found"
    # Return vertex closest to center
    return closest_vertex

def remove_far_to_center(obj: bpy.types.Object, mesh: bmesh.types.BMesh) -> None:
    """Delete long edges and remove unconnected vertices to center cluster"""
    # Keep only vertices connected to the vertex closest to the center
    closest_vertex = get_closest_to_center(mesh.verts)
    # Select closest vertex
    bpy.ops.mesh.select_all(action='DESELECT')
    closest_vertex.select_set(True)
    # Update mesh and free bmesh
    bmesh.update_edit_mesh(obj.data)
    mesh.free()
    # Select linked to vertex and delete other vertices
    bpy.ops.mesh.select_linked()
    bpy.ops.mesh.select_all(action='INVERT')
    bpy.ops.mesh.delete(type='VERT')

def cleanup_section(max_edge_length:float=0.1, method:str="center") -> None:
    """Clean up section.
    Depending on method, remove vertex cluser far from center or non biggest vertex cluster
    """
    obj  = bpy.context.active_object
    mesh = bmesh.from_edit_mesh(obj.data)
    # Remove too long edges
    remove_long_edges(mesh, max_edge_length)
    if method == "center":
        remove_far_to_center(obj, mesh)
    elif method == "cluster":
        # Remove vertices not connected to the biggest edge cluster
        utility.keep_biggest_cluster(mesh, cluster_type="edges")
        # Update mesh and free bmesh
        bmesh.update_edit_mesh(obj.data)
        mesh.free()

def get_section_dimension(obj:bpy.types.Object) -> tuple[np.ndarray, np.ndarray]:
    """Fit a rotating bounding box on set of point and return the dimension of the bounding box"""
    # Set plane as active and switch to Edit mode
    utility.select_all(select=False)
    bpy.context.view_layer.objects.active = obj
    bpy.ops.object.mode_set(mode='EDIT', toggle=False)
    # Siplify geometry to reduce computing time and exclude leaves clusters
#    bpy.ops.mesh.remove_doubles(threshold=0.25)
    cleanup_section(max_edge_length=0.25, method="cluster")
    # Update object data
    obj = bpy.context.active_object
    bpy.ops.object.mode_set(mode='OBJECT', toggle=False)
    bpy.ops.object.mode_set(mode='EDIT', toggle=False)
    # Extract point coordinates and fit bounding box
    vertex_coord = [vertex.co[0:2] for vertex in obj.data.vertices]
    vertex_coord = np.array(vertex_coord)
    print(f"{len(vertex_coord) = }")
    # If only 1 or 0 vertex, cannot fit bounding rectangle, output -1 to indicate error
    if len(vertex_coord) <= 1:
        return tuple(np.array((-1.0,))) * 2
    bbox = minimum_bounding_rectangle(vertex_coord)
    # Add bounding box to active mesh
    draw_polygon_2d(bbox)
    # Measure dimension and center from bounding box coordinates
    try:
        bbox_dim = [np.linalg.norm(bbox[i]-bbox[i+1]) for i in range(2)]
        center = np.average(bbox, axis=0)
    except ValueError:
        # Wrong cross-section, output -1 to indicate error in the cross-section computation
        return tuple(np.array((-1.0,))) * 2
    # Switch back to object mode
    bpy.ops.object.mode_set(mode='OBJECT', toggle=False)
    return (np.array(bbox_dim), np.array(center))

if __name__ == "__main__":
    # Apply cross-section on active object
    test_obj = bpy.context.active_object
    assert test_obj is not None, "No active object"
    test_plane = main(test_obj)
    print(f"{get_section_dimension(test_plane)=}")
