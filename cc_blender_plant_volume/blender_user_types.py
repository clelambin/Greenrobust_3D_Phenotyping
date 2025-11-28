"""List of user defined types used in the blender_plant_volume module"""

# Import libraries
from typing import Literal
from bmesh import types
from mathutils import Vector

# User Literal type
Cartesian = Literal["X", "Y", "Z"]
BooleanOperator = Literal["DIFFERENCE", "INTERSECT"]

# Complex type anotation
BMEntry = types.BMVert | types.BMEdge | types.BMFace
BMFaceList = list[types.BMFace] | None
Segment = tuple[Vector, Vector]
