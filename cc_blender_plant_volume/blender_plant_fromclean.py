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
from collections.abc import Iterable, Callable
from dataclasses import dataclass, field
import numpy as np
import bpy                      # Blender python
import bmesh                    # Blender mesh module
from mathutils import Vector    # Blender object type

# Import user modules
from cc_blender_plant_volume import blender_utility_functions as utility
from cc_blender_plant_volume import blender_extract_attributefiltering as attribute
from cc_blender_plant_volume import blender_point_clustering as cluster
from cc_blender_plant_volume import blender_utility_ransac as ransac
from cc_blender_plant_volume.blender_user_types import Cartesian, BooleanOperator


# Utility functions
def projected_dist(pt1:Vector,
                   pt2:Vector,
                   axis:Cartesian|None="Z") -> float:
    """Compute the distance between 2 points projected on the input axis
    The projection is done by ignoring the given coordinate when computing the distance
    """
    # Set projected axis coordinate to 0
    if axis is not None:
        setattr(pt1, axis.lower(), 0.0)
        setattr(pt2, axis.lower(), 0.0)
    # Compute distance between both points
    return float(np.linalg.norm(np.array(pt1)-np.array(pt2)))

def calc_ratio(nb1:float, nb2:float) -> float:
    """Compute ratio as |nb1 - nb2| / nb1"""
    return abs(nb1-nb2) / nb1


# Classes
@dataclass
class FaceInfo:
    """Geometric characteristic of the current face"""
    area: float
    perimeter: float
    axis_position: float
    center: Vector
    # Field computed in __post_init__
    roundness: float = field(init=False)

    def __post_init__(self) -> None:
        """Compute face roundness
        Roundess computed using the Polsby–Popper applied to a plannar face
        """
        self.roundness = (4 * np.pi * self.area) / (self.perimeter**2)


@dataclass
class FaceCriteria:
    """Set of criterion which define if a face from a section should be considered"""
    min_area: float = float("-inf")
    max_area: float = float("inf")
    min_roundness: float = float("-inf")

    def is_in_range(self, face:FaceInfo) -> bool:
        """Input face match the defined criteria"""
        return (face.area < self.max_area
                and face.area > self.min_area
                and face.roundness > self.min_roundness)


@dataclass
class BranchCriteria:
    """Set of criterion which define if two faces from different section can be connected"""
    max_area_ratio:float = float("inf")
    max_roundness_ratio:float = float("inf")
    max_xy_dist:float = float("inf")
    max_z_dist:float = float("inf")

    def face_in_range(self, source:FaceInfo, target:FaceInfo) -> bool:
        """Check if both faces match criteria"""
        return (calc_ratio(source.roundness, target.roundness) < self.max_roundness_ratio
                and calc_ratio(source.area, target.area) < self.max_area_ratio)

    def distance_in_range(self, xy_dist:float, z_dist:float) -> float|None:
        """Check if horizontal and veritical distance (between faces) in range"""
        return (xy_dist < self.max_xy_dist
                and z_dist < self.max_z_dist)


class FaceNode:
    """Contain information about face"""
    def __init__(self,
                 face:bmesh.types.BMFace,
                 axis_position:float,
                 name:str,
                 save_coords:bool=False) -> None:
        self.name = name
        self.info = FaceInfo(
            area = face.calc_area(),
            perimeter=face.calc_perimeter(),
            center=face.calc_center_median_weighted(),
            axis_position=axis_position
        )
        # If ask to save the cooridnate from the face, extract them from the BMFace
        # Because point coordinate is from a plane initialy normal to Z, the Z coordinate is 0
        # extract only x and y coordinate
        self.coords_2d = [vert.co.xy for vert in face.verts] if save_coords else None
        # At first, parent and children of face node are unknown
        self.parent     = None
        self.children   = set()
        # Mesh vertex id (used to link to created bmesh object)
        self.vertex_id:int|None = None

    def __repr__(self) -> str:
        return f"Face {self.name}"

    def get_closest_face(self,
                         face_list:Iterable["FaceNode"],
                         branch_criteria:BranchCriteria) -> "FaceNode|None":
        """Return the face whose center is closest to current face and with similar roundness"""
        # Initialise closest face index and distance
        closest_dist  = float("inf")
        closest_face: "FaceNode|None" = None
        # Loop through face list and save the face with the closest projected distance
        for face in face_list:
            # Check that faces match ratio criterion
            if not branch_criteria.face_in_range(self.info, face.info):
                continue
            # Compute distance between center and save face if closer than current closest
            xy_dist = projected_dist(self.info.center, face.info.center)
            z_dist  = abs(self.info.axis_position - face.info.axis_position)
            if branch_criteria.distance_in_range(xy_dist, z_dist) and xy_dist < closest_dist:
                closest_dist = xy_dist
                closest_face = face
        return closest_face

    def add_child(self, face:"FaceNode"):
        """Add child to current node if not already in children list"""
        self.children.add(face)

    def add_parent(self, face:"FaceNode"):
        """Add parent to current face and set current face as child of input face"""
        # Can only assign a parent if no parent already assigned
        assert self.parent is None, (
            f"Parent already assigned to {self}\n" \
            f"\tExisting parent: {self.parent}\n"  \
            f"\tAttempt to assign parrent: {face}"
        )
        self.parent = face
        # Add current face as child of parent face
        face.add_child(self)


class TreeStructure:
    """Create tree data structure of face accross multiple sections.

    Using modified vocabulary for the tree data structure:
    - heads : current node without parents
              (once a parent is added to a head, it is removed from this list)
    - leafs : nodes without childres (extremity of tree datastructure)
    - roots : nodes without parents (once all sections have been processed)

    Additionaly, save characteristics for each branches (path from a leaf to a root)
    - branches_xy_dist : sum of horizontal distance between each node of the branch
    """

    def __init__(self,
                 faces:list[FaceNode],
                 branch_criteria:BranchCriteria) -> None:
        """Initialise the tree by defining input nodes as heads and leafs"""
        self.heads:list[FaceNode] = faces
        self.leafs:list[FaceNode] = faces
        self.branch_criteria = branch_criteria

        # Roots are defined after the full tree structure is set
        self.roots:list[FaceNode]|None = None

        # Branch characteristics
        self.branches_xy_dist:list[float] = []

    def __iter__(self):
        """Initialise iterator by setting head to root"""
        if self.roots is None:
            self.set_roots()
        return self.walk_through(self.roots)

    @classmethod
    def walk_through(cls, heads:Iterable[FaceNode]|None):
        """Generator which walk through elements of the tree, from parents to its children"""
        if heads is None:
            return
        for head in heads:
            yield head
            yield from cls.walk_through(head.children)

    @classmethod
    def walk_back(cls, head:FaceNode|None):
        """Generator which walk through elements of the tree, from a child to its parent"""
        if head is None:
            return
        yield head
        yield from cls.walk_back(head.parent)

    def add_section(self, faces:Iterable[FaceNode]) -> None:
        """Look at next section and add faces closest to branches heads"""
        # Initialise the loop
        new_faces:set[FaceNode] = set(faces)
        new_heads:set[FaceNode] = set()
        new_parents:set[FaceNode] = set()
        for head in self.heads:
            closest_face = head.get_closest_face(faces, self.branch_criteria)
            # If no closed faces found for the current head, keep it for the next iteration
            if closest_face is None:
                new_heads.add(head)
                continue
            # If closest face found, set as parent of current face and replace head
            head.add_parent(closest_face)
            new_parents.add(closest_face)
        # Add the faces not added as parents as new leafs
        self.leafs.extend(new_faces.difference(new_parents))
        # Update heads with new head
        self.heads = list(new_heads)
        # Add all faces from the current section to the head
        self.heads.extend(faces)

    def set_roots(self):
        """Set current heads as root"""
        self.roots = self.heads

    def branch_straightness(self, leaf:FaceNode) -> tuple[float, int]:
        """Return the straightness of the branch defined by the input leaf."""
        # Initialise branch info
        nb_nodes  = 0
        branch_xy_dist:float = 0.0
        # Walk to the root, increamenting the horizontal distance
        for node in self.walk_back(leaf):
            # If node has no parent, means it is a root node, so no need to compute distance
            if node.parent is None:
                continue
            node_dist = projected_dist(node.info.center, node.parent.info.center)
            branch_xy_dist += node_dist
            nb_nodes += 1
        return branch_xy_dist, nb_nodes

    def get_straight_branch(self, min_nodes:int=3) -> FaceNode:
        """Walk through the tree structure to compute branch horizontal distance,
        then output leaf of branch with smallest value
        """
        # Initialise loop
        min_dist_ratio  = float("inf")
        min_index = 0
        # Walk through the tree starting from each leaf node to compute branch_xy_dist
        for (index, leaf) in enumerate(self.leafs):
            # Save the branch distance to the class instance
            branch_xy_dist, nb_nodes = self.branch_straightness(leaf)
            self.branches_xy_dist.append(branch_xy_dist)
            # Compare with smalled branch
            if nb_nodes >= min_nodes and branch_xy_dist/nb_nodes < min_dist_ratio:
                min_dist_ratio  = branch_xy_dist/nb_nodes
                min_index = index
        # Get index of branch with minimum horizontal distance
        return self.leafs[min_index]


# Process functions
def remove_ground(obj:bpy.types.Object, ground_height:float=0.002) -> None:
    """Remove all points from mesh below ground height"""
    # Add cube object, shift it toward the ground
    cube_location = Vector([0, 0, ground_height-1])
    bpy.ops.mesh.primitive_cube_add(size=2, location=cube_location)
    ground = utility.get_active_obj()
    ground.name = "Ground"
    # Apply difference boolean to remove anything inside the ground object
    boolean_modifier(source_obj=obj, target_obj=ground, operation="DIFFERENCE", thresh=0)
    # Hide Ground
    utility.hide_object(ground)

def boolean_modifier(source_obj:bpy.types.Object,
                     target_obj:bpy.types.Object,
                     thresh:float=0.00001,
                     operation:BooleanOperator="INTERSECT",
                     apply:bool=True)-> None:
    """Create intersection between object to intersect and plane
    Using boolean mesh operator
    """
    # Mark source object as active
    utility.make_active(source_obj)
    # Boolean modifier using fast intersection mode
    modifier = source_obj.modifiers.new(name="Section", type="BOOLEAN")
    modifier.solver = "FAST"
    modifier.operation = operation
    modifier.double_threshold = thresh
    # Set object to intersect
    modifier.object = target_obj
    # Apply modifier
    if apply:
        bpy.ops.object.modifier_apply(modifier=modifier.name)

def section_faces(obj:bpy.types.Object,
                  axis_position:float,
                  face_criteria:FaceCriteria,
                  min_faces:int=1,
                  save_coords:bool=False) -> list[FaceNode]|None:
    """Remove vertices not connected to any faces from current object mesh,
    Return a list containing the roundness and center of each face in the section
    """
    # Check that input object contain a mesh
    assert isinstance(obj.data, bpy.types.Mesh), f"Object {obj.name} does not contain a mesh"

    # initialise face description list
    face_descr = []

    # Extract mesh from object
    mesh = bmesh.new()
    mesh.from_mesh(obj.data)

    # Save all vertices connected to faces as a set
    connected_verts = set()
    for (index, face) in enumerate(mesh.faces):
        connected_verts.update(face.verts)
        # Save roundness and center of current face
        face_node = FaceNode(face,
                             name = f"{obj.name}-{index}",
                             axis_position = axis_position,
                             save_coords = save_coords)
        # If face roundness and area criterion are respected, add face to the face description
        if face_criteria.is_in_range(face_node.info):
            face_descr.append(face_node)

    # If object has less faces in range than min number of faces, delete it
    if len(face_descr) < min_faces:
        bpy.data.objects.remove(obj, do_unlink=True)
        return None

    # Get vertices not parts of connected vertices and delete them
    all_verts = set(mesh.verts)
    disconnected_verts = list(all_verts.difference(connected_verts))
    utility.delete_vertices(mesh, disconnected_verts)

    # Apply modification to object
    mesh.to_mesh(obj.data)
    mesh.free()

    # Return the face description list, sorted by roundness
    return face_descr

def section_skeleton(section_detail: list[list[FaceNode]],
                     branch_criteria:BranchCriteria) -> TreeStructure:
    """Link face accross sections based on closest projection"""
    # Initialise TreeStructure object containing the link between faces
    section_tree = TreeStructure(section_detail.pop(), branch_criteria)
    # Start from last section, create architecture of linked faces
    for section in reversed(section_detail):
        section_tree.add_section(section)
    return section_tree

def draw_tree(tree_structure:TreeStructure, name:str="PlantStructure") -> None:
    """Draw the tree structure as a bmesh object"""
    # Initialise bmesh object
    mesh = bmesh.new()
    mesh.verts.ensure_lookup_table()
    # Loop through nodes of the tree and creates verts and link to parent node
    for (index, node) in enumerate(tree_structure):
        node_coordinate = node.info.center
        node_coordinate.z = node.info.axis_position
        mesh.verts.new(node_coordinate)
        mesh.verts.ensure_lookup_table()
        node.vertex_id = index
        # If node has parent, create an edge to link both
        if node.parent is not None:
            parent_index = node.parent.vertex_id
            assert parent_index is not None, f"Parent node {node.parent} not added to mesh yet"
            mesh.edges.new((mesh.verts[index], mesh.verts[parent_index]))
    # Create new object and add to current collection
    obj_mesh = bpy.data.meshes.new(name=name)
    obj      = bpy.data.objects.new(name=name, object_data=obj_mesh)
    bpy.context.collection.objects.link(obj)
    # Assign mesh to obj and free mesh
    mesh.to_mesh(obj_mesh)
    mesh.free()

def draw_stick(tree_structure:TreeStructure,
               stick_radius:float=0.003,
               ground_height:float=0.1) -> bpy.types.Object|None:
    """Extract stick from tree structure (straightest branch) and generate a mesh along stick"""
    # Find straightest branch in tree structure
    straight_branch = tree_structure.get_straight_branch()
    # Get coordinates for each node of the branch
    nodes_coord = []
    for node in tree_structure.walk_back(straight_branch):
        # Extract x and y coordinate from the center coordinate and replace z by the section height
        x, y, _ = node.info.center
        nodes_coord.append(Vector([x, y, node.info.axis_position]))
    # If last coordinate higher than ground height, add point using last x, y coordinate
    if node.info.axis_position > ground_height:
        nodes_coord.append(Vector([x, y, ground_height]))

    # Create curve based on nodes_coords
    curve_data = utility.create_curve(nodes_coord, name="Stick")

    # Add bevel (thickness) to the curve
    curve_data.bevel_depth = stick_radius
    curve_data.bevel_mode  = "ROUND"
    curve_data.bevel_resolution = 8
    curve_data.use_fill_caps = True

    # Create object and link to scene
    stick_obj = utility.curve_to_mesh(curve_data)

    # Return created object
    return stick_obj

def reshape_section(sections:dict[Cartesian, FaceNode]) -> list[Vector]:
    """Extract point coordinate from section and reshape the points to the XY quadrant
    Align section to Y and apply symetry to shift points with X<0 to X>0
    """
    # Initialise list containing all points from processed section
    all_points:list[Vector] = []
    # Loop through sections and reshape points if alligned to X
    for (axis, section) in sections.items():
        point_coords = section.coords_2d
        # If coordinate not provided for section, skip it
        if point_coords is None:
            continue
        # If aligned to X, invert the coordinates
        if axis == "X":
            point_coords = [Vector((pt.y, pt.x)) for pt in point_coords]
        # Apply symetry to all points with X<0
        point_coords =  [Vector((abs(pt.x), pt.y)) for pt in point_coords]
        # Append processed points to the list of all points
        all_points.extend(point_coords)
    return all_points

def draw_pot(sections:dict[Cartesian, FaceNode]) -> None:
    """Create a simplified version of the pot, based on the input X and Y sections"""
    # Reshape all sections and collect into a single list of points
    all_points = reshape_section(sections)
    # Use RANSAC to fit a simplified pot section
    ransac_param = ransac.RansacParam(
            nb_sample=5,
            min_cluster=3,
            max_iter=1000,
            dist_thresh=0.0001,
            max_fit=0.9
    )
    init_model = ransac.PotSection()
    fitted_model  = ransac.ransac_fit(init_model, all_points, ransac_param)

    # Print model parameters
    print(f"{fitted_model=}")
    print(f"Ratio of fitted points: {fitted_model.fit_ratio:.3f}")

    # Draw segments
    pot_segments = fitted_model.segments
    for (index, param) in enumerate(fitted_model.params):
        pot_curve    = utility.create_curve(list(pot_segments[index]), name=param.name)
        pot_section  = utility.curve_to_mesh(pot_curve)
        utility.from_z_to_axis(pot_section, "Y")


# TODO: check that obj is updated before getting dimension
def multi_crosssection(obj:bpy.types.Object,
                       face_criteria:FaceCriteria,
                       z_delta:float=0.1) -> list[list[FaceNode]]:
    """Create crosssection of input object spaced by z_delta"""
    # Initialise list which will contain the face information for each section
    section_detail = []
    # Get Z dimention to define number of cross-sections
    obj_dim    = obj.dimensions
    nb_section = int(obj_dim[2] / z_delta) + 1
    for section_index in range(1, nb_section):
        axis_position = section_index * z_delta
        # Create plane at given axis_position height
        plane = utility.create_plane("Z", axis_position)
        # Apply section modifier to plane
        boolean_modifier(source_obj=plane, target_obj=obj)
        # Cleanup cross-section and save list of all faces within section
        face = section_faces(plane, axis_position, face_criteria)
        if face is not None:
            section_detail.append(face)
    # Output face detail for all sections
    return section_detail

def vert_crosssection(obj:bpy.types.Object,
                      face_criteria:FaceCriteria,
                      select:Callable=lambda face:face.info.area) -> dict[Cartesian, FaceNode]:
    """Create vertical crosssection of input object for X and Y direction"""
    # Create a section normal to X and normal to Y
    direction_list:list[Cartesian] = ["X", "Y"]
    face_list = {}
    for direction in direction_list:
        faces = None
        section_dist = 0
        # If the stem intersect with the cross-section, the section fail, retry with offset
        while faces is None and section_dist < 0.02:
            plane = utility.create_plane(direction, section_dist, name_digits=0)
            # Apply section modifier to plane
            boolean_modifier(source_obj=plane, target_obj=obj)
            # Cleanup cross-section and save list of all faces within section
            faces = section_faces(plane, 0, face_criteria, save_coords=True)
            section_dist += 0.001
        # Add face matching the input select function to the dictionary
        if faces is not None:
            face_list[direction] = sorted(faces, key=select)[-1]
    # Output face detail
    return face_list


def plant_cleanup(plant:bpy.types.Object) -> None:
    """Extract green plant from 3D model, remove support and ouput model attribute
    Model attribute: plant volume, plant height, area projection
    """
    # Set criterion requirement for face and branch
#    plant_criteria = FaceCriteria(
#           min_area = 1e-6,
#           max_area = 1e-3,
#           min_roundness = 0.5
#    )
    pot_criteria = FaceCriteria(
            min_area = 1e-4,
            max_area = 1e-2
    )
    stick_criteria = FaceCriteria(
            min_area = 5e-6,
            max_area = 3e-5,
            min_roundness = 0.5
    )
    branch_criteria = BranchCriteria(
            max_area_ratio = 1.0,
            max_roundness_ratio = 1.0,
            max_xy_dist = 0.005,
            max_z_dist = 0.1
    )
    # Check if any vertices below ground
    min_z = min(obj_bbox[2] for obj_bbox in plant.bound_box)
    if min_z < 0:
        # Remove vertices below ground
        remove_ground(plant)

    # Create vertical section, used to extract the pot
    pot_section = vert_crosssection(plant, pot_criteria)
    print(f"{pot_section=}")
    draw_pot(pot_section)
#    # Use attribute filtering to exclude pot and get green plant
#    plant_green = attribute.exclude_pot(plant)
#    deleted_vertices = cluster.dbscan_filter(plant_green)
#    # If -1 returned as number of deleted vertices, no cluster detected, raise an error
#    assert deleted_vertices != -1, f"No cluster detected for {plant.name}"
#
#    # Find stick (straightest branch) and draw cordesponding cylinder
#    draw_stick(tree_structure)


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
