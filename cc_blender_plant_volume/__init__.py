"""Blender volume computation for digitised 3D plant model.
Can be called on a signle model (by defining MODEL_FILE) or on a whole folder.
Position and align plant model to the origin, scale the plant based on the plant pot dimension,
then compute the volume of the green part of the plant.

Workflow:
- Import obj
- Remove background noise
- Use color attribute to detect pot and cup
- Position pot center to global center
- Allign object center to +Z axis
- Intersect obj with XY plane to get pot cross-section
- Use pot cross-section to scale model
- Compute volume of green plant
"""

# Libraries
import os             # File manager
import bpy            # Blender python
# User module
from cc_blender_plant_volume import blender_plant_alignment as modelprep

# User variables
WORKING_PATH   = r"C:\Users\cleme\Documents\Hohenheim\00_Courses\320_Landscape_and_Plant_Ecology\MSc_3D_Plant_Characterisation\3D_Digitalisation\02_Input_Greenhouse\2025-05-13_Harvest_Nicotina_benthamiana"
OUTPUT_PATH    = r"C:\Users\cleme\Documents\Hohenheim\00_Courses\320_Landscape_and_Plant_Ecology\MSc_3D_Plant_Characterisation\3D_Digitalisation\04_Output_BlenderScript"
OUTPUT_TABLE   = "Plant_data.csv"
MODEL_FOLDER   = "02_Metashape_BatchProc_Conf"
MODEL_FILE_EXT = "ply"
MODEL_FILE     = None
#MODEL_FILE     = "Metashape_NB004_S30_20250612_2048.ply"
#MODEL_FILE     = "Metashape_NB005_S34_20250612_2121.ply"

# Operator classes
# (Registering operators to be callable in Blender)
class PlantModelPrep(bpy.types.Operator):
    """Operator for single plant model preparation.
    Allign plant model, apply scaling and compute volume.
    """
    # Blender identifiers
    bl_idname  = "object.plant_model_prep"
    bl_label   = "Model preparation for single plant"
    bl_options = {"REGISTER", "UNDO"}

    # Class properties
    model_file_ext: bpy.props.EnumProperty(
        name = "File extension type",
        description = "File extension type",
        items = (
            # name, value, description
            ("ply", "ply", "Standford PLY (.ply)"),
            ("obj", "obj", "Wavefront (.obj)"),
        ),
        default = "ply",
    )
    working_path  : bpy.props.StringProperty(default=".")
    model_folder  : bpy.props.StringProperty()
    model_file    : bpy.props.StringProperty()

    def pool(self, context):
        """Check if context suitable to call the function"""
        return bpy.context.object.mode == "OBJECT"

    def execute(self, context):
        """Execute operator"""
        model_path = os.path.join(self.working_path, self.model_folder, self.model_file)
        modelprep.single_model_prep(model_path, file_ext=self.model_file_ext)
        return {"FINISHED"}

class PlantFolderLookup(bpy.types.Operator):
    """Operator running plant model preparation for all folders in input path.
    Allign plant model, apply scaling and compute volume.
    """
    # Blender identifiers
    bl_idname  = "object.plant_folder_lookup"
    bl_label   = "Model preparation for whole folder"
    bl_options = {"REGISTER", "UNDO"}

    # Class properties
    model_file_ext: bpy.props.EnumProperty(
        name = "File extension type",
        description = "File extension type",
        items = (
            # name, value, description
            ("ply", "ply", "Standford PLY (.ply)"),
            ("obj", "obj", "Wavefront (.obj)"),
        ),
        default = "ply",
    )
    working_path  : bpy.props.StringProperty(default=".")
    model_folder  : bpy.props.StringProperty()
    output_path   : bpy.props.StringProperty()
    output_table  : bpy.props.StringProperty()

    def pool(self, context):
        """Check if context suitable to call the function"""
        return bpy.context.object.mode == "OBJECT"

    def execute(self, context):
        """Execute operator"""
        modelprep.loop_through_folders(self.working_path,
                                       self.model_folder,
                                       self.model_file_ext,
                                       self.output_path,
                                       self.output_table)
        return {"FINISHED"}


# Register full classes
classes = (
    PlantModelPrep,
    PlantFolderLookup,
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

    # If MODEL_FILE is None, loop through whole folder, otherwise, process only specified file
    if MODEL_FILE is None:
        bpy.ops.object.plant_folder_lookup(
            working_path   = WORKING_PATH,
            model_folder   = MODEL_FOLDER,
            model_file_ext = MODEL_FILE_EXT,
            output_path    = OUTPUT_PATH,
            output_table   = OUTPUT_TABLE,
        )
    else:
        bpy.ops.object.plant_model_prep(
            working_path   = WORKING_PATH,
            model_folder   = MODEL_FOLDER,
            model_file     = MODEL_FILE,
            model_file_ext = MODEL_FILE_EXT,
        )
