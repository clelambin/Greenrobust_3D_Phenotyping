"""Main script for the model processing from already aligned and scaled model

Note:
- require digitisation with marker for alignment
- replace model_prep from blender_plant_modelprep.py

Workflow:
- Import obj
- Create vertical section to detect pot shape and horizontal section to detect stick shape
- Create simplified pot shape based on the section and remove it from the plant model
- Create cylinder at location of wood stick and remove it from the plant model
- Extract model attribute (volume, ...)
"""

# Import libraries
import operator
from collections.abc import Iterable, Callable
from dataclasses import dataclass, field
from typing import Literal
import os
import numpy as np              # Array manipulation
import bpy                      # Blender python
import bmesh                    # Blender mesh module
from mathutils import Matrix, Vector    # Blender object type

# Import user modules
from cc_blender_plant_volume import blender_utility_functions as utility
from cc_blender_plant_volume import blender_point_clustering as cluster
from cc_blender_plant_volume import blender_utility_ransac as ransac
from cc_blender_plant_volume import blender_plant_metrics as metrics
from cc_blender_plant_volume import blender_plant_rendering as render
from cc_blender_plant_volume.blender_user_types import Cartesian, ConfigOptions

# Script function
FILE_PREFIX = "Metashape_"
SPECIES_CONFIG: dict[str, dict[ConfigOptions, float|bool]] = {
    # Individual pot
    "AT":{"POT_OFFSET": 0, "HAS_STICK":False, "POT_WIDTH": 0.075, "POT_HEIGHT": 0.092, "POT_THICK": 0.013},
    "BD":{"POT_OFFSET": 0, "HAS_STICK":True,  "POT_WIDTH": 0.075, "POT_HEIGHT": 0.092, "POT_THICK": 0.013},
    "BR":{"POT_OFFSET": 0, "HAS_STICK":True,  "POT_WIDTH": 0.095, "POT_HEIGHT": 0.11, "POT_THICK": 0.013},
    "HS":{"POT_OFFSET": 0, "HAS_STICK":True,  "POT_WIDTH": 0.095, "POT_HEIGHT": 0.11, "POT_THICK": 0.013},
    "HV":{"POT_OFFSET": 0, "HAS_STICK":True,  "POT_WIDTH": 0.095, "POT_HEIGHT": 0.11, "POT_THICK": 0.013},
    "NB":{"POT_OFFSET": 0, "HAS_STICK":False, "POT_WIDTH": 0.12,  "POT_HEIGHT": 0.128, "POT_THICK": 0.013},
    "SD":{"POT_OFFSET": 0, "HAS_STICK":True,  "POT_WIDTH": 0.138, "POT_HEIGHT": 0.138, "POT_THICK": 0.013},
    "SL":{"POT_OFFSET": 0, "HAS_STICK":False, "POT_WIDTH": 0.138, "POT_HEIGHT": 0.142, "POT_THICK": 0.013},
    "TA":{"POT_OFFSET": 0, "HAS_STICK":True,  "POT_WIDTH": 0.095, "POT_HEIGHT": 0.11, "POT_THICK": 0.013},
    # Population pots
    "ATp":{"POT_OFFSET": 0, "HAS_STICK":False, "POT_WIDTH": 0.21, "POT_HEIGHT": 0.10, "POT_THICK": 0.05},
    "BDp":{"POT_OFFSET": 0, "HAS_STICK":False, "POT_WIDTH": 0.21, "POT_HEIGHT": 0.10, "POT_THICK": 0.05},
    "BRp":{"POT_OFFSET": 0, "HAS_STICK":False, "POT_WIDTH": 0.30, "POT_HEIGHT": 0.13, "POT_THICK": 0.05},
    "HSp":{"POT_OFFSET": 0, "HAS_STICK":False, "POT_WIDTH": 0.30, "POT_HEIGHT": 0.13, "POT_THICK": 0.05},
    "HVp":{"POT_OFFSET": 0, "HAS_STICK":False, "POT_WIDTH": 0.30, "POT_HEIGHT": 0.13, "POT_THICK": 0.05},
    "NBp":{"POT_OFFSET": 0, "HAS_STICK":False, "POT_WIDTH": 0.36, "POT_HEIGHT": 0.13, "POT_THICK": 0.05},
    "SDp":{"POT_OFFSET": 0, "HAS_STICK":False, "POT_WIDTH": 0.41, "POT_HEIGHT": 0.15, "POT_THICK": 0.05},
    "SLp":{"POT_OFFSET": 0, "HAS_STICK":False, "POT_WIDTH": 0.41, "POT_HEIGHT": 0.15, "POT_THICK": 0.05},
    "TAp":{"POT_OFFSET": 0, "HAS_STICK":False, "POT_WIDTH": 0.30, "POT_HEIGHT": 0.13, "POT_THICK": 0.05},
}
SPECIES_OVERWRITE: dict[str, dict[ConfigOptions, float|bool]] = {
    # Square pots (higher pot height)
    "BR03IC": {"POT_HEIGHT": 0.123, "POT_THICK": 0.015},
    "BR07IC": {"POT_HEIGHT": 0.123, "POT_THICK": 0.015},
    "BR08IC": {"POT_HEIGHT": 0.123, "POT_THICK": 0.015},
    "BR16ID": {"POT_HEIGHT": 0.123, "POT_THICK": 0.015},
    "BR17ID": {"POT_HEIGHT": 0.123, "POT_THICK": 0.015},
    "BR18ID": {"POT_HEIGHT": 0.123, "POT_THICK": 0.015},
    "BR21IC": {"POT_HEIGHT": 0.123, "POT_THICK": 0.015},
    "BR22IC": {"POT_HEIGHT": 0.123, "POT_THICK": 0.015},
    "BR23IC": {"POT_HEIGHT": 0.123, "POT_THICK": 0.015},
    "NB23IC": {"POT_HEIGHT": 0.152, "POT_THICK": 0.015},
    "NB26IC": {"POT_HEIGHT": 0.152, "POT_THICK": 0.015},
    "NB29IC": {"POT_HEIGHT": 0.152, "POT_THICK": 0.015},
    "SD09IC": {"POT_HEIGHT": 0.152, "POT_THICK": 0.015},
}
CONFIG_DEFAULT:dict[ConfigOptions, float|bool] = {"POT_OFFSET": 0.003, "HAS_STICK":True}


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
        using the Polsby–Popper applied to a planar face
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
        """Check if horizontal and vertical distance (between faces) in range"""
        return (xy_dist < self.max_xy_dist
                and z_dist < self.max_z_dist)

@dataclass
class ModifierParam:
    """Boolean Operator setting to apply via boolean_modifier function"""
    operator:Literal["DIFFERENCE", "INTERSECT"] = "INTERSECT"
    vertex_criteria:Literal["RANGE", "BBOX", "NONE"] = "NONE"
    vertex_ratio_range:tuple[float, float]|None=None
    thresh:float=0.00001
    apply:bool=True


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
        # If ask to save the coordinate from the face, extract them from the BMFace
        # Because point coordinate is from a plane initially normal to Z, the Z coordinate is 0
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
        """Return the face whose centre is closest to current face and with similar roundness"""
        # Initialise closest face index and distance
        closest_dist  = float("inf")
        closest_face: "FaceNode|None" = None
        # Loop through face list and save the face with the closest projected distance
        for face in face_list:
            # Check that faces match ratio criterion
            if not branch_criteria.face_in_range(self.info, face.info):
                continue
            # Compute distance between centre and save face if closer than current closest
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
    """Create tree data structure of face across multiple sections.

    Using modified vocabulary for the tree data structure:
    - heads : current node without parents
              (once a parent is added to a head, it is removed from this list)
    - leafs : nodes without childres (extremity of tree datastructure)
    - roots : nodes without parents (once all sections have been processed)

    Additionally, save characteristics for each branches (path from a leaf to a root)
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
        # Walk to the root, incrementing the horizontal distance
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
            # Compare with smaller branch
            if nb_nodes >= min_nodes and branch_xy_dist/nb_nodes < min_dist_ratio:
                min_dist_ratio  = branch_xy_dist/nb_nodes
                min_index = index
        # Get index of branch with minimum horizontal distance
        return self.leafs[min_index]


# Process functions
def remove_ground(obj:bpy.types.Object, ground_height:float=0.005) -> None:
    """Remove all points from mesh below ground height"""
    # Add cube object, shift it toward the ground
    cube_location = Vector([0, 0, ground_height-1])
    bpy.ops.mesh.primitive_cube_add(size=2, location=cube_location)
    ground = utility.get_active_obj()
    ground.name = "Ground"
    # Apply difference boolean to remove anything inside the ground object
    boolean_modifier(source_obj=obj,
                     target_obj=ground,
                     param=ModifierParam(operator="DIFFERENCE", thresh=0))
    # Hide Ground
    utility.hide_object(ground)

def modified_vertex_count(obj:bpy.types.Object) -> int:
    """Return the number of vertex after modifier are applied to input object"""
    depsgraph = bpy.context.evaluated_depsgraph_get()
    evaluated_obj = obj.evaluated_get(depsgraph)
    assert isinstance(evaluated_obj.data, bpy.types.Mesh)
    return len(evaluated_obj.data.vertices)

def is_in_range(value:float, range_tuple:tuple[float, float]) -> bool:
    """Return true if the value is within the input range"""
    return range_tuple[0] < value < range_tuple[1]

def boolean_modifier(source_obj:bpy.types.Object,
                     target_obj:bpy.types.Object,
                     param:ModifierParam=ModifierParam()) -> float:
    """Create intersection between object to intersect and plane
    Using boolean mesh operator, return vertex ratio (modified/initial vertex count)
    """
    # Get vertex count before modification
    assert isinstance(source_obj.data, bpy.types.Mesh)
    vertex_init = len(source_obj.data.vertices)
    # Mark source object as active
    utility.make_active(source_obj)
    # Boolean modifier using fast intersection mode
    modifier = source_obj.modifiers.new(name="Section", type="BOOLEAN")
    assert isinstance(modifier, bpy.types.BooleanModifier)
    modifier.solver = "FAST"
    modifier.operation = param.operator
    modifier.double_threshold = param.thresh
    # Set object to intersect
    modifier.object = target_obj

    # If using auto vertex criteria, define vertex ratio range
    # based number of source vertices in target object bounding box
    if param.vertex_criteria == "BBOX":
        # Count vertices from source in bounding box of target
        source_verts = source_obj.data.vertices
        vertex_in_bbox = len(utility.get_verts_in_bbox(source_verts, target_obj))
        # Update vertex ratio range based on vertex count in target object and operator type
        if param.operator == "INTERSECT":
            # In intersection, expect target vertices to be around
            # number of vertices in target object (+/- 20%)
            param.vertex_ratio_range = (
                    vertex_in_bbox / vertex_init * 0.8,
                    vertex_in_bbox / vertex_init * 1.2
                    )
        elif param.operator == "DIFFERENCE":
            # In difference, expect target vertices to be all but the ones in target object
            # (in general, expect more remaining vertices than present in bounding box
            # since bounding box bigger than actual object)
            param.vertex_ratio_range = (
                    (vertex_init - vertex_in_bbox) / vertex_init * 0.9,
                    (vertex_init - vertex_in_bbox) / vertex_init * 1.35
                    )

    # If vertex_ratio_range specified, only apply modifier if within range
    if param.vertex_ratio_range is not None:
        # Evaluate the vertex count after modifier modification
        vertex_out = modified_vertex_count(source_obj)
        if not is_in_range(vertex_out/vertex_init, param.vertex_ratio_range):
            # Switch to Exact solver and allow self intersection
            modifier.solver = "EXACT"
            modifier.use_self = True
        # If the vertex range is still not in range, return -1 and disable modification
        vertex_out = modified_vertex_count(source_obj)
        if not is_in_range(vertex_out/vertex_init, param.vertex_ratio_range):
            modifier.show_viewport = False
            modifier.show_render = False
            return -1

    # Apply modifier
    if param.apply:
        bpy.ops.object.modifier_apply(modifier=modifier.name)
    # Return vertex ratio
    return modified_vertex_count(source_obj) / vertex_init

def section_faces(obj:bpy.types.Object,
                  axis_position:float,
                  face_criteria:FaceCriteria,
                  min_faces:int=1,
                  save_coords:bool=False) -> list[FaceNode]|None:
    """Remove vertices not connected to any faces from current object mesh,
    Return a list containing the roundness and centre of each face in the section
    """
    # Check that input object contain a mesh
    assert isinstance(obj.data, bpy.types.Mesh), f"Object {obj.name} does not contain a mesh"

    # initialise face description list
    face_descr = []

    # Save all vertices connected to faces as a set
    with utility.bmesh_edit(obj) as mesh:
        connected_verts = set()
        for (index, face) in enumerate(mesh.faces):
            connected_verts.update(face.verts)
            # Save roundness and centre of current face
            face_node = FaceNode(face,
                                 name = f"{obj.name}-{index}",
                                 axis_position = axis_position,
                                 save_coords = save_coords)
            # If face roundness and area criterion are respected, add face to the face description
            if face_criteria.is_in_range(face_node.info):
                face_descr.append(face_node)
        # Get vertices not parts of connected vertices and delete them
        all_verts = set(mesh.verts)
        disconnected_verts = list(all_verts.difference(connected_verts))
        utility.delete_vertices(mesh, disconnected_verts)

    # If object has less faces in range than min number of faces, delete it
    if len(face_descr) < min_faces:
        bpy.data.objects.remove(obj, do_unlink=True)
        return None

    # Return the face description list, sorted by roundness
    return face_descr

def section_skeleton(section_detail: list[list[FaceNode]],
                     branch_criteria:BranchCriteria) -> TreeStructure|None:
    """Link face across sections based on closest projection"""
    # If section_detail is empty, return None to prevent script to stop
    if len(section_detail) == 0:
        return None
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
               ground_height:float=0.1,
               hide:bool=True) -> bpy.types.Object|None:
    """Extract stick from tree structure (straightest branch) and generate a mesh along stick"""
    # Find straightest branch in tree structure
    straight_branch = tree_structure.get_straight_branch()
    # Get coordinates for each node of the branch
    nodes_coord = []
    for node in tree_structure.walk_back(straight_branch):
        # Extract x and y coordinate from the centre coordinate and replace z by the section height
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

    # If requested, hide stick
    if hide:
        utility.hide_object(stick_obj)

    # Return created object
    return stick_obj

def reshape_section(sections:dict[Cartesian, FaceNode]) -> list[Vector]:
    """Extract point coordinate from section and reshape the points to the XY quadrant
    Align section to Y and apply symmetry to shift points with X<0 to X>0
    """
    # Initialise list containing all points from processed section
    all_points:list[Vector] = []
    # Loop through sections and reshape points if aligned to X
    for (axis, section) in sections.items():
        point_coords = section.coords_2d
        # If coordinate not provided for section, skip it
        if point_coords is None:
            continue
        # If aligned to X, invert the coordinates
        if axis == "X":
            point_coords = [Vector((pt.y, pt.x)) for pt in point_coords]
        # Apply symmetry to all points with X<0
        point_coords =  [Vector((abs(pt.x), pt.y)) for pt in point_coords]
        # Append processed points to the list of all points
        all_points.extend(point_coords)
    return all_points

# TODO: not very clean way to extract parameter from PotSection
def create_pot(pot_section:ransac.PotSection,
               offset:float=0.005,
               top_tickness:float=0.001,
               name:str="Simplified_pot") -> bpy.types.Object:
    """Create simplified pot model based on fitted pot section"""
    # Extract all parameters from fitted model and convert to dictionary
    param_dict = pot_section.get_param()
    # Extract relevant parameters as variables
    pot_width_base = param_dict["pot_width"][0]
    pot_width_top = param_dict["pot_width"][1]
    pot_height = param_dict["pot_height"][0]
    soil_width = param_dict["soil_width"][1]
    soil_height = param_dict["soil_height"][0]

    # Start from a primitive cube
    bpy.ops.mesh.primitive_cube_add()
    # Save primitive cube as object
    cube = bpy.context.active_object
    assert cube is not None, "Could not create primitive cube"
    cube.name = name

    # Edit the primitive to fit parameters
    with utility.bmesh_edit(cube) as cube_mesh:
        # Select lower face and rescale to fit base and top width
        base_verts = [vert for vert in cube_mesh.verts if vert.co.z == -1]
        top_verts = [vert for vert in cube_mesh.verts if vert.co.z == 1]
        bmesh.ops.transform(cube_mesh, matrix=Matrix.Scale(pot_width_base, 3), verts=base_verts)
        bmesh.ops.transform(cube_mesh, matrix=Matrix.Scale(pot_width_top, 3), verts=top_verts)
        for verts in base_verts:
            verts.co.z = 0
        for verts in top_verts:
            verts.co.z = pot_height
        # Inset top face to translate downward to fit soil height
        top_face = utility.extract_face_from_verts(top_verts)
        assert top_face is not None, "No faces connected to input vertices"
        # If specify top_thickness, add initial inset to create a thickness at the pot top
        # The top face remain the middle one after inset,
        # so can use the same face over multiple operations
        if top_tickness>0:
            bmesh.ops.inset_individual(cube_mesh, faces=[top_face], thickness=top_tickness)
        # Second inset used to create the soil face (translated downward)
        bmesh.ops.inset_individual(cube_mesh, faces=[top_face], thickness=pot_width_top-soil_width)
        bmesh.ops.translate(cube_mesh,
                            vec=Vector((0, 0, soil_height-pot_height)),
                            verts=list(top_face.verts))

        # Offset all faces
        utility.offset_faces(cube_mesh, dist=offset)
    # Return modified cube
    return cube

def fit_pot(sections:dict[Cartesian, FaceNode],
            hide_pot:bool=True,
            pot_offset:float=0.005) -> tuple[bpy.types.Object, ransac.PotSection]:
    """Fit simplified pot, based on the input X and Y sections,
    Return created object and fitted pot model
    """
    # Reshape all sections and collect into a single list of points
    all_points = reshape_section(sections)
    # Use RANSAC to fit a simplified pot section
    ransac_param = ransac.RansacParam(
            nb_sample=5,
            max_iter=10000,
            min_pts_per_line=0,
            dist_thresh=0.001,
            max_fit=0.6
    )
    init_model = ransac.PotSection()
    fitted_model  = ransac.ransac_fit(init_model, all_points, ransac_param)

    # Print model parameters
    print(f"{fitted_model=}")
    print(f"Ratio of fitted points: {fitted_model.fit_ratio:.3f}")

#    # Draw segments
#    pot_segments = fitted_model.segments
#    for (index, param) in enumerate(fitted_model.params):
#        pot_curve    = utility.create_curve(list(pot_segments[index]), name=param.name)
#        pot_section  = utility.curve_to_mesh(pot_curve)
#        utility.from_z_to_axis(pot_section, "Y")
#
#    point_cluster = fitted_model.cluster_points(all_points, ransac_param.dist_thresh)
#    print(f"Points per segments: {[len(cluster) for cluster in point_cluster]}")
#    print(f"Segments: {pot_segments}")

    # Draw create pot object based on fitted pot section
    pot = create_pot(fitted_model, offset=pot_offset)

    if hide_pot:
        utility.hide_object(pot)
    return pot, fitted_model

def multi_crosssection(obj:bpy.types.Object,
                       face_criteria:FaceCriteria,
                       z_delta:float=0.1,
                       hide_plane:bool=True) -> list[list[FaceNode]]:
    """Create cross-section of input object spaced by z_delta"""
    # Initialise list which will contain the face information for each section
    section_detail = []
    # Get the highest point in Z to define number of cross-sections
    assert isinstance(obj.data, bpy.types.Mesh), f"Object {obj.name} does not have a mesh"
    max_z = max(vertex.co[2] for vertex in obj.data.vertices)
    nb_section = int(max_z / z_delta) + 1
    for section_index in range(1, nb_section):
        axis_position = section_index * z_delta
        # Create plane at given axis_position height
        plane = utility.create_plane("Z", axis_position, hide=hide_plane)
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
                      select:Callable=lambda face:face.info.area,
                      hide_plane:bool=True) -> dict[Cartesian, FaceNode]:
    """Create vertical cross-section of input object for X and Y direction"""
    # Create a section normal to X and normal to Y
    direction_list:list[Cartesian] = ["X", "Y"]
    face_list = {}
    for direction in direction_list:
        faces = None
        section_dist = 0
        # If the stem intersect with the cross-section, the section fail, retry with offset
        while faces is None and section_dist < 0.02:
            plane = utility.create_plane(direction, section_dist, name_digits=0, hide=hide_plane)
            # Apply section modifier to plane
            boolean_modifier(source_obj=plane, target_obj=obj, param=ModifierParam(thresh=0))
            # Cleanup cross-section and save list of all faces within section
            faces = section_faces(plane, 0, face_criteria, save_coords=True)
            section_dist += 0.001
        # Add face matching the input select function to the dictionary
        if faces is not None:
            face_list[direction] = sorted(faces, key=select)[-1]
    # Output face detail
    return face_list


def plant_cleanup(plant:bpy.types.Object,
                  pot_offset:float=0.005,
                  remove_stick:bool=True,
                  pot_thickness:float=0.001,
                  population:bool=False,
                  pot_dim:None|tuple[float, float]=None) -> tuple[Vector, int]:
    """Extract green plant from 3D model, remove support and output model attribute"""
    # Initialise code error
    # If no error, code error remain at 0, otherwise, each error adds up to it
    # (similar to chmod ownership code)
    # code error and corresponding function
    # +1 : pot boolean operation
    # +2 : stick boolean operation
    # +4 : pot RANSAC shape did not reach requirement
    code_error:int = 0

    # Set criterion requirement for face and branch
#    plant_criteria = FaceCriteria(
#           min_area = 1e-6,
#           max_area = 1e-3,
#           min_roundness = 0.5
#    )
    pot_criteria = FaceCriteria(
            min_area = 1e-4,
            max_area = 1e-1
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

    # If not input pot dimension, fit pot
    if pot_dim is None:
        # Create vertical section, used to extract the pot
        pot_section = vert_crosssection(plant, pot_criteria)
    #    print(f"{pot_section=}")
    #    for (axis, section) in pot_section.items():
    #        section_curve = utility.create_curve(section.coords_2d, name=f"Section_{str(axis)}")
    #        section_mesh  = utility.curve_to_mesh(section_curve)
    #        utility.from_z_to_axis(section_mesh, axis)
        pot, pot_param = fit_pot(pot_section, pot_offset=pot_offset)

        # If did not fit defined fitted ratio, update code error
        if not pot_param.reach_fit:
            code_error += 4
    else:
        # Use input pot dimensions
        pot_param = ransac.preset_pot(pot_dim[0], pot_dim[1])
        pot = create_pot(pot_param, offset=pot_offset, top_tickness=pot_thickness)

    # Define plant reference point at Z=soil_height
    plant_ref = Vector((0, 0, pot_param.get_param()["soil_height"][0]))

    # Remove pot from plant
    pot_vert_ratio = boolean_modifier(plant,
                                      pot,
                                      param=ModifierParam(
                                          operator="DIFFERENCE",
                                          vertex_criteria="BBOX"))
    # If vertex ratio is -1, boolean not applied, update code error
    if pot_vert_ratio == -1:
        code_error += 1

    if remove_stick:
        # Detecting the wood stick from the plant by creating horizontal sections
        # and looking at section structure
        section_detail = multi_crosssection(plant, stick_criteria, z_delta=0.01)
        tree_structure = section_skeleton(section_detail, branch_criteria)
#       draw_tree(tree_structure)

        # Find stick (straightest branch) and draw corresponding cylinder
        stick = draw_stick(tree_structure) if tree_structure is not None else None
        # Remove stick from the plant (if stick detected)
        if stick is not None:
            stick_vert_ratio = boolean_modifier(plant,
                                                stick,
                                                param=ModifierParam(
                                                    operator="DIFFERENCE",
                                                    vertex_ratio_range=(0.3, 1)))
            # If vertex ratio is -1, boolean not applied, update code error
            if stick_vert_ratio == -1:
                code_error += 2

    # Run DBScan clustering on vertex to remove vertex cluster further from given distance
    # Number of cluster depend on population vs individual pot
    nb_cluster = 9 if population else 1
    deleted_vertices = cluster.dbscan_filter(plant, dbscan_eps=0.02, nb_cluster=nb_cluster)
    # If -1 returned as number of deleted vertices, no cluster detected, raise an error
    assert deleted_vertices != -1, f"Not enough cluster detected for {plant.name}"

    # Return plant reference points and code error (used for metrics)
    return plant_ref, code_error

def plant_metrics(plant:bpy.types.Object,
                  output_dir:str|None=None,
                  plant_ref:Vector|None=None) -> dict[str, float]:
    """Extract metrics from input plant model and save rendered images for quality control
    Model attribute: plant volume, plant height, area projection
    """
    # If reference point specified, convert it to NumPy array
    plant_ref_array = np.array(plant_ref) if plant_ref is not None else None
    # Compute plant metrics (volume, surface and dimensions)
    # (If output dir specified, add temp image in output dir for area project)
    tmp_img = os.path.join(output_dir, "ObjectProject.png") if output_dir is not None else None
    metrics_dict = metrics.calc_metrics(plant, tmp_img, ref=plant_ref_array)

    # Return computed metrics
    print(f"{metrics_dict = }")

    # If output_dir specified, save rendered image of prepared model in output directory
    if output_dir is not None:
        render.model_rendering(output_dir, plant.name)

    # Cleanup unused data
    bpy.ops.outliner.orphans_purge()

    # Return computed metrics
    return metrics_dict

def single_model_prep(obj_path:str|None=None,
                      output_path:str|None=None,
                      file_ext:str="obj") -> dict[str, float]:
    """Single model preparation: import the model, compute the volume and save output blend file"""

    # Prepare working environment
    # If input obj_path, import it, otherwise, work on selected object
    if obj_path is not None:
        utility.cleanup_env()
        model = utility.import_file(obj_path, file_ext)
        # Read species label from path
        file_name = os.path.basename(obj_path)
    else:
        model = bpy.context.active_object
        # If no active object alert the user
        assert model is not None, "No active object"
        # Read species from object name
        file_name = model.name

    # Extract relevant species config from model name
    species_short = file_name.replace(FILE_PREFIX, "")[0:2]

    # Check if pot is individual or population pot
    species_orga = file_name.replace(FILE_PREFIX, "")[3:4]
    # If population pot, append p at end of species short to indicate it
    if species_orga == "P":
        species_short = species_short + "p"

    # Find configuration for corresponding species
    species_config = SPECIES_CONFIG.get(species_short, CONFIG_DEFAULT)

    # If plant listed in overwrite, overwrite defined values
    plant_label = file_name.replace(FILE_PREFIX, "")[0:6]
    plant_overwrite = SPECIES_OVERWRITE.get(plant_label, {})
    for param, value in plant_overwrite.items():
        species_config[param] = value

    # Prepare the model and compute all plant metrics
    # Extract pot dimensions if found
    pot_dim = (species_config.get("POT_WIDTH", None), species_config.get("POT_HEIGHT", None))
    pot_dim = None if any(dim is None for dim in pot_dim) else pot_dim

    # Cleanup plant model
    plant_ref, code_error = plant_cleanup(model,
                                          species_config["POT_OFFSET"],
                                          species_config["HAS_STICK"],
                                          pot_thickness=species_config["POT_THICK"],
                                          pot_dim=pot_dim,
                                          population=species_orga=="P")

    # Read plant metrics
    model_metrics = plant_metrics(model, output_dir=output_path, plant_ref=plant_ref)
    # Update metrics code error to output from plant_cleanup function
    model_metrics["Code_Error"] = code_error
    # Save blend file in output folder
    if output_path is not None:
        utility.save_blend(obj_path, output_path)
    # Return plant metrics
    return model_metrics

def process_model_folder(file_list:list[str],
                         working_path:str,
                         volume_path:str,
                         output_path:str="",
                         file_ext:str="obj",
                         metrics_keys:list[str]=["Volume"]) -> None:
    """Loop through file list and process all model files"""
    for name in file_list:
        if not name.lower().endswith(file_ext) or "ptscloud" in name.lower():
            continue
        # Process plant model
        model_path = os.path.join(working_path, name)
        model_metrics = single_model_prep(model_path, output_path=output_path, file_ext=file_ext)
        # Save all plant metrics to csv file
        plant_info = [working_path, name]
        for key in metrics_keys:
            plant_info.append(str(model_metrics.get(key, 0)))
        with open(volume_path, "a", encoding="utf-8") as volume_file:
            volume_file.write(f"{','.join(plant_info)}\n")

def loop_through_folders(working_path:str,
                         model_folder:str,
                         model_file_ext:str,
                         output_path:str,
                         output_table:str) -> None:
    """Loop through the folders and process folders matching input name"""
    # Loop through all files and process plant model
    volume_path = os.path.join(output_path, output_table)
    plant_header = ["Path", "Plant name"]
    metrics_keys = list(metrics.calc_metrics().keys())
    with open(volume_path, "w", encoding="utf-8") as volume_file:
        volume_file.write(f"{','.join(plant_header)},{','.join(metrics_keys)}\n")
    for root, _, files in os.walk(working_path):
        # Check if last folder in root matches with specified model folder
        last_folder = root.split(os.sep)[-1]
        if model_folder == last_folder:
            print(f"Processing {root}")
            process_model_folder(files,
                                 root,
                                 volume_path,
                                 output_path,
                                 model_file_ext,
                                 metrics_keys=metrics_keys)

def main() -> None:
    """Main function, prepare model for selected model on the 3D view"""
    model_metrics = single_model_prep()
    print(f"{model_metrics=}")

if __name__ == "__main__":
    # Test model preparation on active object
    main()
