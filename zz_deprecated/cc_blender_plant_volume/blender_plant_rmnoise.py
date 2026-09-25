"""Reduce model noise from plant model.
Expend mesh and run self-intersection to check for nearby meshes, remove the isolated elements
(Use attribute and vertex index to compare processed model with initial, lead to unreliable results)
"""

import bpy
import numpy as np

# To do:
# - complete type anotation (use bpy.types instead of list[float])

# Warning:
# - bpy.ops tend to slow down scripts, to check for alternative when possible

# Blender mesh operations
def reduce_mesh_count(dist:float=0.01) -> None:
    """Use select double to remove nearby vertices to reduce mesh count of selected object"""
    # Switch to edit mode
    bpy.ops.object.mode_set(mode='EDIT', toggle=False)
    # Decimate all faces
    bpy.ops.mesh.select_all(action='SELECT')
    # Decimate does not preserve attribute, use remove double instead
#    bpy.ops.mesh.decimate(ratio=ratio)
    bpy.ops.mesh.remove_doubles(threshold=dist)

def expand_mesh(dist:float=0.1) -> None:
    """Expend mesh of selected object(used to connect element"""
    # Switch to edit mode
    bpy.ops.object.mode_set(mode='EDIT', toggle=False)
    # Apply shrink flatten on all faces
    bpy.ops.mesh.select_all(action='SELECT')
    bpy.ops.transform.shrink_fatten(value=dist)

def self_intersect() -> None:
    """Self-intersect faces of selected object"""
    # Switch to edit mode
    bpy.ops.object.mode_set(mode='EDIT', toggle=False)
    # Self-intersect all faces
    bpy.ops.mesh.select_all(action='SELECT')
    bpy.ops.mesh.intersect(mode='SELECT', separate_mode='NONE', solver='FAST')
    # Unselect intersection
    bpy.ops.mesh.select_all(action="DESELECT")

def distance_to_center(point_coord:list[float], center_coord:list[float]) -> float:
    """Compute the distance from a given point to a center coordinate"""
    # Convert coordinates to numpy array
    point_array = np.array(point_coord)
    center_array = np.array(center_coord)
    # Return distance between both points
    return np.linalg.norm(point_array - center_array)

def first_point_in_range(point_list:list[list[float]], center:list[float], threshold:float) -> int:
    """Return the index of the first point within a threshold distance from the center. If no point in range, return the closest instead"""
    # Initialise the loop
    # (float("inf") initialise the distance to infinite, so the first distance will always be smaller than that)
    closest_to_center = -1
    min_distance = float("inf")
    # Loop through point and stop if reach a point within the distance threshold
    for index, point in enumerate(point_list):
        distance = distance_to_center(point, center)
        if distance < threshold:
            # Distance within threshold, return the index and stop the loop
            return index
        elif distance < min_distance:
            # Distance smaller than current min distance, update the closest point
            min_distance = distance
            closest_to_center = index
    # Return the closest point to center
    print("No dstance within threshold, return closest distance instead")
    print(f"Threshold distance: {threshold}, min distance: {min_distance}")
    return closest_to_center

def reset_center() -> None:
    """Reset object location to its current center location (needed after moving)"""
    # Only work in Object mode
    bpy.ops.object.mode_set(mode='OBJECT', toggle=False)
    # Set origin to the selected object
    bpy.ops.object.origin_set(type = 'ORIGIN_GEOMETRY')

# Not really a list of float
def local_to_global(matrix, coord) -> list[float]:
    """Convert vertex coordinate from local to global"""
    # With numpy @ act as a matrix multiplier
    return matrix @ coord

def select_from_center(center_range_ratio:float=0.05) -> None:
    """Select a vertex in the center of the object bounding box and select vertices linked to this one"""
    # Switch to edit mode
    bpy.ops.object.mode_set(mode='EDIT', toggle=False)
    # Extract object center, dimension and list of vertices coordinates
    obj = bpy.context.active_object
    obj_center = obj.location
    obj_size   = obj.dimensions
    # Original point coord are in local coordinate, convert to global
#    obj_vertices = [vertex.co for vertex in obj.data.vertices]
    obj_matrix = obj.matrix_world
    obj_vertices = [local_to_global(obj_matrix, point.co) for point in obj.data.vertices]
    # Scale the range distance to the center based on object dimension
    range_dist = center_range_ratio * np.linalg.norm(np.array(obj_size))
    # Extract index of a point close to the center
    near_center = first_point_in_range(obj_vertices, obj_center, range_dist)
    # Debug
    print(f"Selected index: {near_center}")
    print(f"Corresponding coordinate: {obj_vertices[near_center]}")
    print(f"Center coordinate: {obj_center}")
    # Select corresponding point and expend selection to select linked vertices to that point
    select_from_index(obj, [near_center])
#    # (To select a point from list, nead to be in object mode)
#    bpy.ops.object.mode_set(mode = 'OBJECT', toggle=False)
#    obj.data.vertices[near_center].select = True
#    bpy.ops.object.mode_set(mode='EDIT', toggle=False)
#    bpy.ops.mesh.select_linked()

def select_from_index(obj, vertex_list:list[int], select_link:bool=True) -> None:
    """Select vertices from index list and linked vertices"""
    # Unselect all
    bpy.ops.object.mode_set(mode = 'EDIT', toggle=False)
    bpy.ops.mesh.select_all(action="DESELECT")
    # (To select a point from list, nead to be in object mode)
    bpy.ops.object.mode_set(mode = 'OBJECT', toggle=False)
    # Create numpy array which is set to True for each indices in the list
    vertex_nb = len(obj.data.vertices)
    to_select = np.zeros(vertex_nb, dtype=bool)
    to_select[np.array(vertex_list)] = True
    # Select vertices from boolean array and connected
    obj.data.vertices.foreach_set("select", to_select)
    bpy.ops.object.mode_set(mode='EDIT', toggle=False)
    if select_link:
        bpy.ops.mesh.select_linked()

# list output can be of different type based on attribute
def get_attribute_from_selected(attribute_name:str) -> list:
    """Return a list of attributes from all selected vertices"""
    # Inspired from https://blender.stackexchange.com/questions/1412/efficient-way-to-get-selected-vertices-via-python-without-iterating-over-the-en
    # Reading attribute only work in Object mode
    bpy.ops.object.mode_set(mode='OBJECT', toggle=False)
    obj = bpy.context.active_object
    vertex_nb = len(obj.data.vertices)
    # Check attribute type
    attribute_type = obj.data.attributes[attribute_name].data_type
    # Initialise which sellected vertices and relevant attribute for each vertices
    selected_vertices = np.zeros(vertex_nb, dtype=bool)
    attribute_vertices = np.zeros(vertex_nb, dtype=attribute_type.lower())
    # For each vertices, get if selected and relevant attribute
    obj.data.vertices.foreach_get('select', selected_vertices)
    obj.data.attributes[attribute_name].data.foreach_get("value", attribute_vertices)
    # Ignore all 0 values in the attribute: all points added later are initialise with a 0 value for the attribute
    attribute_nonzero = attribute_vertices != 0
    # Return non zero attributes for selected vertices
    return attribute_vertices[selected_vertices & attribute_nonzero]

# select class
# (not sure if class relevant here)
class select_mesh:
    """Create a copy of an object with offset and self-intesect to select the plant and related components"""
    def __init__(self, obj : bpy.types.Object) -> None:
        """Read the object and save mesh info"""
        # Save active object
        self.obj = obj
        self.obj_data = bpy.context.collection.objects[obj.name].data

    def duplicate_object(self) -> None:
        """Duplicate the object and save original vertex id as attributes of duplicated object"""
        # Duplicate object
        bpy.ops.object.duplicate()
        self.obj_dup = bpy.context.active_object
        self.obj_dup_data = bpy.context.collection.objects[self.obj_dup.name].data
        # Save vertex id as attribute
        dup_attr = self.obj_dup_data.attributes.new(name="source_id", type="INT", domain="POINT")
        # (because vertex not edited, vertex id is 0:len(vertices)
        vertex_id = range(0, len(self.obj_data.vertices))
        dup_attr.data.foreach_set("value", vertex_id)

def main() -> None:
    """Main function, remove background noise on selected object"""
    # Initialise class and working object
    reset_center()
    obj = bpy.context.active_object
    select_proc = select_mesh(obj)

    # Duplicate object and work on duplicated object
    select_proc.duplicate_object()
    reduce_mesh_count()
    expand_mesh()
    self_intersect()
    select_from_center()
    source_vertices = get_attribute_from_selected(attribute_name="source_id")
#    for vertices in source_vertices:
#        print(vertices)

    # Work back on original object and select vertices from duplicated object source id
    bpy.context.view_layer.objects.active = obj
    # Remove dupplicated and select original object
    select_proc.obj_dup.select_set(False)
    bpy.data.objects.remove(select_proc.obj_dup)
    select_proc.obj.select_set(True)
    select_from_index(obj, source_vertices, select_link=True)
    # Inverse selection and delete selected
    bpy.ops.object.mode_set(mode='EDIT', toggle=False)
    bpy.ops.mesh.select_all(action='INVERT')
    bpy.ops.mesh.delete(type='VERT')

if __name__ == "__main__":
    main()
