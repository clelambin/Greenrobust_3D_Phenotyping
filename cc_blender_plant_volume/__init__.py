"""Blender volume computation for digitised 3D plant model.
Can be called on a single model (by defining MODEL_FILE) or on a whole folder.
Position and align plant model to the origin, scale the plant based on the plant pot dimension,
then compute the volume of the green part of the plant.

Workflow:
- Import obj.
- Remove background noise.
- Use colour attribute to detect pot and cup.
- Position pot centre to global centre.
- Align object centre to +Z axis.
- Intersect obj with XY plane to get pot cross-section.
- Use pot cross-section to scale model.
- Compute volume of green plant.
"""

# Libraries
import os             # File manager
import bpy            # Blender python
# User module
from cc_blender_plant_volume import blender_plant_fromclean as modelprep
from cc_blender_plant_volume import blender_utility_functions as utility

# Operator classes
# (Registering operators to be callable in Blender)
class PlantModelPrep(bpy.types.Operator):
    """Operator for plant model preparation.
    Align plant model, apply scaling and compute volume.
    If no MODEL_FILE specified in input config, loop through all folders in input path.
    """
    # Blender identifiers
    bl_idname  = "object.plant_model_prep"
    bl_label   = "Plant model preparation"
    bl_options = {"REGISTER", "UNDO"}

    def pool(self, context):
        """Check if context suitable to call the function"""
        return bpy.context.object.mode == "OBJECT"

    def execute(self, context):
        """Execute operator"""
        modelprep.process_from_config()
        return {"FINISHED"}


# Register full classes
classes = (
    PlantModelPrep,
)

def register():
    """Register specified classes to Blender"""
    for _class in classes:
        bpy.utils.register_class(_class)

def unregister():
    """Unregister specified classes to Blender"""
    for _class in classes:
        bpy.utils.unregister_class(_class)

if __name__ == "__main__":
    register()

    # Load config file and run model preparation
    bpy.ops.object.plant_model_prep()
