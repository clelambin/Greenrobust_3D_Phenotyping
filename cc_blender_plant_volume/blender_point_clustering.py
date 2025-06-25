# Import libraries
import open3d as o3d  # Used for point clustering
import numpy as np    # Array and matrix operations
import bpy            # Blender python
import bmesh          # Blender mesh module

# Import user modules
from cc_blender_plant_volume import blender_utility_functions as utility

def dbscan_clustering(vertices_coord, dbscan_eps=0.1) -> np.ndarray:
    """Use dbscan to cluster vertices list
    Return label indicating cluster number as int, -1 represent isolated vertices
    """
    # Convert vertices numpy array to point cloud
    vertices_ptcloud = o3d.geometry.PointCloud()
    vertices_ptcloud.points = o3d.utility.Vector3dVector(vertices_coord)
    with o3d.utility.VerbosityContextManager(o3d.utility.VerbosityLevel.Debug) as _:
        labels = np.array(vertices_ptcloud.cluster_dbscan(eps=dbscan_eps,
                                                          min_points=10,
                                                          print_progress=True))
    # Return labels indicating the cluster value
    return labels

def assign_attribute(obj:bpy.types.Object, labels:np.ndarray) -> None:
    """Assign label to input object as new attribute"""
    # Switch to object mode to edit point attribute
    bpy.ops.object.mode_set(mode='OBJECT', toggle=False) # Go to object mode
    # Add label as vertices atribute
    vertices_attribute = obj.data.attributes.new(name="dbscan_label", type="INT", domain="POINT")
    vertices_attribute.data.foreach_set("value", labels)

def get_biggest_cluster(label:np.ndarray) -> np.intp:
    """Return integer of the cluster containing the most amount of points"""
    # Remove all -1 values from label (represent isolated point)
    label_noneg = np.delete(label, np.where(label == -1)[0])
    # If no negative value found, return -1 and alert the user
    if len(label_noneg) == 0:
        print("Warning, empty sequence, only isolated points found")
        return np.intp(-1)
    # Count the frequency of each remaining cluster label
    label_count = np.bincount(label_noneg)
    # Return most frequent value
    return np.argmax(label_count)

def exclude_label(mesh:bmesh.types.BMesh,
                  target_label:np.intp) -> list[bmesh.types.BMVert]:
    """Return list of vertices not including given label"""
    # Load vertex list and saved attribute
    vertex_list = mesh.verts
    attribute   = mesh.verts.layers.int["dbscan_label"]
    # Return vertex whose attribute does not match the target
    return [vertex for vertex in vertex_list if vertex[attribute] != target_label]

def dbscan_filter(obj:bpy.types.Object) -> int:
    """Apply dbscan clustering on input object and keep biggest cluster.
    Return number of deleted vertices
    """
    bpy.ops.object.mode_set(mode='OBJECT', toggle=False) # Go to object mode
    # Read object vertex coordinates and convert to np_array
    vertices_coord = np.array([vertex.co for vertex in obj.data.vertices])
    # If no vertices found, empty object, return -1 to indicate that no cluster found
    if len(vertices_coord) == 0: return -1
    # Run dbscan clustering
    cluster_label = dbscan_clustering(vertices_coord)
    print(f"{len(cluster_label)=}")
    assign_attribute(obj, cluster_label)
    # Get label of biggest cluster
    biggest_cluster = get_biggest_cluster(cluster_label)
    # If no cluster found, return -1 to indicate that no cluster found
    if biggest_cluster == -1: return -1
    # Switch to edit mesh and load bmesh
    utility.make_active(obj)
    bpy.ops.object.mode_set(mode='EDIT', toggle=False)
    obj_mesh = bmesh.from_edit_mesh(obj.data)
    # Delete vertices not belonging to the biggest cluster
    vertex_toexclude = exclude_label(obj_mesh, biggest_cluster)
    utility.delete_vertices(obj_mesh, vertex_toexclude)
    # Switch back to object mode
    bpy.ops.object.mode_set(mode='OBJECT', toggle=False)
    # Return number of deleted vertices
    return len(vertex_toexclude)

def ransac_plane(obj: bpy.types.Object, ransac_dist=0.1) -> np.ndarray:
    """Return coordinate of the plane fitting point cluster (using RANSAC algorythm)"""
    bpy.ops.object.mode_set(mode='OBJECT', toggle=False) # Go to object mode
    # Read object vertex coordinates and convert to np_array
    vertices_coord = np.array([vertex.co for vertex in obj.data.vertices])
    # Load vertice corrdinate as point cloud
    vertices_ptcloud = o3d.geometry.PointCloud()
    vertices_ptcloud.points = o3d.utility.Vector3dVector(vertices_coord)
    # Fit plane to point cloud
    plane_model, _ = vertices_ptcloud.segment_plane(distance_threshold=ransac_dist,
                                                    ransac_n=3,
                                                    num_iterations=1000)
    # Return plane coordinate
    return plane_model

if __name__ == "__main__":
    # Test dbscan filtering on active object
    active_obj = bpy.context.active_object
    assert active_obj is not None, "No active object"
    dbscan_filter(active_obj)
