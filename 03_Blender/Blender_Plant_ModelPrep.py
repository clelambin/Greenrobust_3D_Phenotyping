# ModelPrep: main script for the 3d model processing

# User variables
module_path  = r"C:\Users\cleme\Documents\Hohenheim\00_Courses\320_Landscape_and_Plant_Ecology\MSc_3D_Plant_Characterisation\3D_Digitalisation\00_Scripts\Blender"
working_path = r"C:\Users\cleme\Documents\Hohenheim\00_Courses\320_Landscape_and_Plant_Ecology\MSc_3D_Plant_Characterisation\3D_Digitalisation\02_Input_Greenhouse\2025-05-13_Harvest_Nicotina_benthamiana"
output_path  = r"C:\Users\cleme\Documents\Hohenheim\00_Courses\320_Landscape_and_Plant_Ecology\MSc_3D_Plant_Characterisation\3D_Digitalisation\04_Output_BlenderScript"
volume_file  = "Plant_volume.csv"
#working_path = r"C:\Users\cleme\Documents\Hohenheim\00_Courses\320_Landscape_and_Plant_Ecology\MSc_3D_Plant_Characterisation\3D_Digitalisation\02_Input_Greenhouse\2025-05-13_Harvest_Nicotina_benthamiana\01_Metashape_BatchProc"
#model_file   = "Metashape_NB005_S26_20250519_1428.obj"

# Workflow:
# - Import obj
# - Remove background noise
# - Generate plant skeleton to extract pot center
# - Position pot center to global center
# - Allign object center to +Z axis
# - Intersect obj with XY plane to get pot cross-section
# - Use pot cross-section to scale model

# Import libraries
import bpy            # Blender python
import bmesh          # Blender mesh module
import mathutils      # Blender object type
import os             # File manager
import sys            # Edit python path
import numpy as np    # Array and math operations

# Edit python path to import user modules
sys.path.insert(0, module_path)

# Import user modules
import Blender_Extract_Skeleton
import Blender_Plant_RmNoise
import Blender_Extract_CrossSection

def cleanup_env(obj_to_remove:list[str] = ["Cube"], type_to_remove:list[str] = ["MESH"]):
    """Remove non-relevant object from scene"""
    # Unselect all to avoid deleting previously selected object
    bpy.ops.object.select_all(action='DESELECT')
    
    # Loop through object and delete object in the list
    for obj in bpy.data.objects:
        if obj.type in type_to_remove or obj.name in obj_to_remove:
            obj.select_set(True)
            bpy.ops.object.delete(use_global=False)

def import_obj(filepath:str):
    """Import obj file"""
    # Check if file exist and is obj
    if not os.path.isfile(filepath):
        raise OSError(f"File {filepath} not found")
    if not filepath.endswith(".obj"):
        raise OSError(f"File {filepath} is not an obj")
        
    # Import file
    bpy.ops.wm.obj_import(filepath=filepath)
    
def get_location(obj=None) -> mathutils.Vector:
    """Reset center and get location of active object or sepified object"""
    bpy.ops.object.origin_set(type='ORIGIN_GEOMETRY')
    if obj is None:
        obj = bpy.context.active_object
    return obj.location

def translate_to_center(ref_obj:bpy.types.Object) -> None:
    """Move all mesh objects to center reference object"""
    bpy.ops.object.mode_set(mode='OBJECT', toggle=False)
    # Set reference object as active
    bpy.ops.object.select_all(action='DESELECT')
    ref_obj.select_set(True)
    # Get reference object location
    ref_location = get_location()
    # Translate all mesh object based on ref location
    bpy.ops.object.select_all(action='DESELECT')
    for obj in bpy.data.objects:
        if obj.type == 'MESH':
            obj.select_set(True)
            bpy.ops.transform.translate(value=-ref_location, orient_type='GLOBAL')
            obj.select_set(False)
            
def allign_to_z(ref_obj:bpy.types.Object) -> None:
    """Rotate object around 2 axis to allign center of ref object to Z axis"""
    bpy.ops.object.mode_set(mode='OBJECT', toggle=False)
    # Mark reference object as active
    bpy.ops.object.select_all(action='DESELECT')
    bpy.context.view_layer.objects.active = ref_obj
    # Set pivot point to center of grid
    bpy.context.scene.cursor.location = mathutils.Vector((0, 0, 0))
    bpy.context.scene.tool_settings.transform_pivot_point = 'CURSOR'
    # Perform the two rotation
    for obj in bpy.data.objects:
        if obj.type == 'MESH':
            obj.select_set(True)
            # - Rotation along Y
            center = get_location(ref_obj)
            rotation_angle = np.pi/2 - np.arctan2(center.z, center.x)
            bpy.ops.transform.rotate(value=rotation_angle,
                                     orient_axis='Y',
                                     orient_type='GLOBAL',
                                     center_override=(0.0, 0.0, 0.0))
            # - Rotation along X
            center = get_location()
            rotation_angle = -np.pi/2 + np.arctan2(center.z, center.y)
            bpy.ops.transform.rotate(value=rotation_angle,
                                     orient_axis='X',
                                     orient_type='GLOBAL',
                                     center_override=(0.0, 0.0, 0.0))
            obj.select_set(False)
            
def calc_scale_ratio(measure:np.ndarray, expected:float, method:str) -> float:
    """Return scale ratio to scale to the expected value"""
    avg_measure = getattr(np, method, "average")(measure)
    return expected / avg_measure
            
def scale_to_ref(ratio) -> None:
    """Scale all mesh object to defined ratio"""
    # Initialise transformation
    bpy.ops.object.mode_set(mode='OBJECT', toggle=False)
    bpy.ops.object.select_all(action='DESELECT')
    bpy.context.scene.tool_settings.transform_pivot_point = 'CURSOR'
    # Loop through all object and scale mesh objects
    for obj in bpy.data.objects:
        if obj.type == 'MESH':
            obj.select_set(True)
            bpy.ops.transform.resize(value=(ratio, ratio, ratio), orient_type='GLOBAL',
                                     center_override=(0.0, 0.0, 0.0))
            obj.select_set(False)
    
# Warning: if switch scaling to bmesh, would need to update        
def calc_volume(obj:bpy.types.Object) -> float:
    """Return volume from bmesh after application of transformation"""
    # Mark object as active and switch to Edit mode
    # ERROR: Operator bpy.ops.object.select_all.poll() failed, context is incorrect
    #        (Happened in Edit mode)
    bpy.ops.object.mode_set(mode='OBJECT', toggle=False)
    bpy.ops.object.select_all(action='DESELECT')
    bpy.context.view_layer.objects.active = obj
    bpy.ops.object.mode_set(mode='EDIT', toggle=False)
    # Extract bmesh from object
    orig_mesh = bmesh.from_edit_mesh(obj.data)
    temp_mesh = orig_mesh.copy()
    # Apply mesh transformation
    temp_mesh.transform(obj.matrix_world)
    # Calculate volume
    volume = temp_mesh.calc_volume()
    # Free up bmesh and switch back to Object mode
    orig_mesh.free()
    temp_mesh.free()
    bpy.ops.object.mode_set(mode='OBJECT', toggle=False)
    # Return volume
    return volume

def import_model(obj_path:str) -> None:
    """Prepare the working environment and load the 3D model"""
    # Prepare working environment
    cleanup_env()
    import_obj(obj_path)
    
def save_blend(obj_path:str, output_path:str) -> None:
    """Save the prepared model as blend file in output file"""
    # Extract plant name
    plant_name = obj_path.split(os.sep)[-1]
    plant_name = plant_name.split(".")[0]
    plant_file = f"{plant_name}.blend"
    bpy.ops.wm.save_as_mainfile(filepath=os.path.join(output_path, plant_file))

def model_prep(pot_size:float=0.13) -> float:
    """Main function, prepare the geometry and output the prepared volume"""
    plant = bpy.context.active_object
    
    # Clean up obj
    Blender_Plant_RmNoise.main()
    Blender_Extract_Skeleton.main(name="Pot")
    Blender_Extract_Skeleton.delete_isolated()
    Blender_Extract_Skeleton.keep_biggest_cluster()
    
    # Allign Object based on pot center
    pot = bpy.context.active_object
    translate_to_center(ref_obj = pot)
    allign_to_z(plant)
    
    # Compute Cross-section
    section = Blender_Extract_CrossSection.main(plant)
    section_dim = Blender_Extract_CrossSection.get_section_dimension(section)
    
    # Perform scaling based on measured ratio if the cross-section returned proper values
    if len(section_dim) != 1:
        ratio = calc_scale_ratio(section_dim, pot_size, method="min")
        scale_to_ref(ratio)
        volume_coef = 1
    else:
        # section returned -1, skip the scaling and multiply the volume by -1 to alert the user
        volume_coef = -1
    
    # Compute volume
    plant_volume = calc_volume(plant) * volume_coef
    print(f"{plant_volume = }")
    
    # Cleanup unused data
    bpy.ops.outliner.orphans_purge()
    
    # Return the volume
    return plant_volume

def single_model_prep(obj_path:str, output_path:str="", pot_size:float=0.13) -> float:
    """Single model preparation: import the model, compute the volume and save output blend file"""
    # Import the model
    import_model(obj_path)
    # Prepare the model and compute the volume
    volume = model_prep(pot_size)
    # Save blend file in output folder
    if output_path != "":
        save_blend(obj_path, output_path)
    # Return the volume
    return volume

def loop_through_files(file_list:list[str], volume_path:str) -> None:
    """Loop through file list and process all obj files"""
    for name in files:
        if not name.lower().endswith(".obj") or "ptscloud" in name.lower():
            continue
        # Process plant model
        model_path = os.path.join(root, name)
        volume = single_model_prep(model_path, output_path=output_path)
        # Save volume to csv file
        with open(volume_path, "a") as volume_file:
            volume_file.write(f"{root},{name},{volume}\n")
    
if __name__ == "__main__":
#    # Test on one file
#    model_path = os.path.join(working_path, model_file)
#    volume = single_model_prep(model_path, output_path=output_path)
    
    # Loop through all files and process plant model
    volume_path = os.path.join(output_path, volume_file)
    with open(volume_path, "w") as volume_file:
        volume_file.write("Path,Plant name,Volume\n")
    for root, dirs, files in os.walk(working_path):
        if not "01_Metashape_BatchProc" in root:
            continue
        print(f"Processing {root}")
        loop_through_files(files, volume_path)