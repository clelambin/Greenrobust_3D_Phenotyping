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
from typing import Literal
from collections.abc import Iterable
from dataclasses import dataclass, field
import numpy as np
import bpy            # Blender python
import bmesh          # Blender mesh module
import mathutils      # Blender object type

# Import user modules
from cc_blender_plant_volume import blender_utility_functions as utility

# User type
Cartesian = Literal["X", "Y", "Z"]
BooleanOperator = Literal["DIFFERENCE", "INTERSECT"]


# Utility functions
def projected_dist(pt1:mathutils.Vector,
                   pt2:mathutils.Vector,
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
    section_z: float
    center: mathutils.Vector
    # Field computed in __post_init__
    roundness: float = field(init=False)

    def __post_init__(self) -> None:
        """Compute face roundness
        Roundess computed using the Polsby–Popper applied to a plannar face
        """
        self.roundness = (4 * np.pi * self.area) / (self.perimeter**2)

class FaceNode:
    """Contain information about face"""
    def __init__(self,
                 face:bmesh.types.BMFace,
                 section_z:float|None = None,
                 name:str|None=None) -> None:
        self.name = name
        self.info = FaceInfo(
            area = face.calc_area(),
            perimeter=face.calc_perimeter(),
            center=face.calc_center_median(),
            section_z=section_z
        )
        # At first, parent and children of face node are unknown
        self.parent     = None
        self.children   = set()
        # Mesh vertex id (used to link to created bmesh object)
        self.vertex_id:int|None = None

    def __repr__(self) -> str:
        return f"Face {self.name}"

    def get_closest_face(self,
                         face_list:Iterable["FaceNode"],
                         area_ratio:float = 0.9,
                         roundness_ratio:float = 0.9) -> "FaceNode|None":
        """Return the face whose center is closest to current face and with similar roundness"""
        # Initialise closest face index and distance
        closest_dist  = float("inf")
        closest_face: "FaceNode|None" = None
        # Loop through face list and save the face with the closest projected distance
        for face in face_list:
            # Skip faces with roundness or area too far from current face
            if calc_ratio(self.info.roundness, face.info.roundness) > roundness_ratio:
                continue
            if calc_ratio(self.info.area, face.info.area) > area_ratio:
                continue
            # Compute distance between center and save face if closer than current closest
            dist = projected_dist(self.info.center, face.info.center, axis="Z")
            if dist < closest_dist:
                closest_dist = dist
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

    def __init__(self, faces:list[FaceNode]) -> None:
        """Initialise the tree by defining input nodes as heads and leafs"""
        self.heads:list[FaceNode] = faces
        self.leafs:list[FaceNode] = faces

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
            closest_face = head.get_closest_face(faces)
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
            try:
                node_dist = projected_dist(node.info.center, node.parent.info.center)
                branch_xy_dist += node_dist
                nb_nodes += 1
            except AttributeError:
                pass
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
    cube_location = mathutils.Vector([0, 0, ground_height-1])
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

def create_plane(normal_axis:Cartesian="Z",
                 axis_position:float=0,
                 name_digits:int=3) -> bpy.types.Object:
    """Create a plane normal to given axis and passing by normal axis as input position"""
    # Define unit vetor used to define position for each axis
    unit_vect = {
        "X": mathutils.Vector([1,0,0]),
        "Y": mathutils.Vector([0,1,0]),
        "Z": mathutils.Vector([0,0,1]),
    }
    # Define rotation matrix for each axis
    rotation_mat = {
        "X": mathutils.Matrix([[0,0,1],[0,1,0],[-1,0,0]]),
        "Y": mathutils.Matrix([[1,0,0],[0,0,-1],[0,1,0]]),
        "Z": mathutils.Matrix([[1,0,0],[0,1,0],[0,0,1]]),
    }
    # Create plane at given location
    # (by default the plane is alligned with z)
    location = axis_position*unit_vect[normal_axis]
    bpy.ops.mesh.primitive_plane_add(size=10, align='WORLD', location=location)
    plane = utility.get_active_obj()
    # .{name_digit}f add digit to the number to improve object order
    plane.name = f"Plane_{normal_axis}_{axis_position:.{name_digits}f}"
    # Return error if plane could not be created
    assert plane is not None, "Error, section plane not created"
    # Rotate plane to be normal to specified axis
    plane.rotation_euler.rotate(rotation_mat[normal_axis])
    # Return plane object
    return plane

def section_faces(obj:bpy.types.Object,
                  section_z:float,
                  min_faces:int=2,
                  min_roundness:float=0.5,
                  min_area:float=1e-7) -> list[FaceNode]|None:
    """Remove vertices not connected to any faces from current object mesh,
    Return a list containing the roundness and center of each face in the section
    """
    # initialise face description list
    face_descr = []

    # Extract mesh from object
    mesh = bmesh.new()
    mesh.from_mesh(obj.data)

    # If object has less faces than min number of faces, delete it (do not process)
    if len(mesh.faces) < min_faces:
        bpy.data.objects.remove(obj, do_unlink=True)
        return None

    # Save all vertices connected to faces as a set
    connected_verts = set()
    for (index, face) in enumerate(mesh.faces):
        connected_verts.update(face.verts)
        # Save roundness and center of current face
        face_node = FaceNode(face, name=f"{obj.name}.{index}", section_z=section_z)
        # If face roundness and are above threshold, add it to the face description
        if face_node.info.roundness > min_roundness and face_node.info.area > min_area:
            face_descr.append(face_node)

    # Get vertices not parts of connected vertices and delete them
    all_verts = set(mesh.verts)
    disconnected_verts = list(all_verts.difference(connected_verts))
    utility.delete_vertices(mesh, disconnected_verts)

    # Apply modification to object
    mesh.to_mesh(obj.data)
    mesh.free()

    # Return the face description list, sorted by roundness
    return sorted(face_descr, key=lambda face: face.info.roundness, reverse=True)

def section_skeleton(section_detail: list[list[FaceNode]]) -> TreeStructure:
    """Link face accross sections based on closest projection"""
    # Initialise TreeStructure object containing the link between faces
    section_tree = TreeStructure(section_detail.pop())
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
        node_coordinate.z = node.info.section_z
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
    mesh.to_mesh(obj.data)
    mesh.free()

def draw_stick(tree_structure:TreeStructure, stick_radius:float=0.003) -> bpy.types.Object:
    """Extract stick from tree structure (straightest branch) and generate a mesh along stick"""
    # Find straightest branch in tree structure
    straight_branch = tree_structure.get_straight_branch()
    # Get coordinates for each node of the branch
    nodes_coord = []
    for node in tree_structure.walk_back(straight_branch):
        # Extract x and y coordinate from the center coordinate and replace z by the section height
        x, y, _ = node.info.center
        nodes_coord.append(mathutils.Vector([x, y, node.info.section_z]))

    # Initialise the curve object
    curve_data = bpy.data.curves.new("Stick", "CURVE")
    curve_data.dimensions = "3D"
    # Add polyline to curve with a point at center of each node of straight branch
    curve_poly = curve_data.splines.new("POLY")
    curve_poly.points.add(len(nodes_coord) - 1)
    for (index, coord) in enumerate(nodes_coord):
        # Polyline coodinate is composed of 4 components: x, y, z (part of coord) and the weight
        curve_poly.points[index].co = (*coord, 1)

    # Add bevel (thickness) to the curve
    curve_data.bevel_depth = stick_radius
    curve_data.bevel_mode  = "ROUND"
    curve_data.bevel_resolution = 8
    curve_data.use_fill_caps = True

    # Create object and link to scene
    curve_obj  = bpy.data.objects.new("Stick", curve_data)
    curve_mesh = bpy.data.meshes.new_from_object(curve_obj)
    stick_obj  = bpy.data.objects.new("Stick_meshed", curve_mesh)
    bpy.context.collection.objects.link(curve_obj)
    bpy.context.collection.objects.link(stick_obj)

    # Return created object
    return stick_radius

# TODO: check that obj is updated before getting dimension
def multi_crosssection(obj:bpy.types.Object, z_delta:float=0.1) -> list:
    """Create crosssection of input object spaced by z_delta"""
    # Initialise list which will contain the face information for each section
    section_detail = []
    # Get Z dimention to define number of cross-sections
    obj_dim    = obj.dimensions
    nb_section = int(obj_dim[2] / z_delta) + 1
    for section_index in range(1, nb_section):
        section_z = section_index * z_delta
        # Create plane at given section_z height
        plane = create_plane("Z", section_z)
        # Apply section modifier to plane
        boolean_modifier(source_obj=plane, target_obj=obj)
        # Cleanup cross-section and save list of all faces within section
        face = section_faces(plane, section_z)
        if face is not None:
            section_detail.append(face)
    # Output face detail for all sections
    return section_detail

def plant_cleanup(obj:bpy.types.Object) -> None:
    """Extract green plant from 3D model, remove support and ouput model attribute
    Model attribute: plant volume, plant height, area projection
    """
    # Check if any vertices below ground
    min_z = min(obj_bbox[2] for obj_bbox in obj.bound_box)
    if min_z < 0:
        # Remove vertices below ground
        remove_ground(obj)

    # Plant section and extract face detail for each section
    section_detail = multi_crosssection(obj, z_delta=0.01)
    tree_structure   = section_skeleton(section_detail)
    draw_tree(tree_structure)

    # Find stick (straightest branch) and draw cordesponding cylinder
    draw_stick(tree_structure)

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
