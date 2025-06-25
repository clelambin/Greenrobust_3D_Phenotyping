# ModelPrep: main script for the 3d model processing

# Workflow:
# - Import obj
# - Remove background noise
# - Generate plant skeleton to extract pot center
# - Position pot center to global center
# - Allign object center to +Z axis
# - Intersect obj with XY plane to get pot cross-section
# - Use pot cross-section to scale model

# Import libraries
import os             # File manager
import bpy            # Blender python
import bmesh          # Blender mesh module
import mathutils      # Blender object type
import numpy as np    # Array and matrix operations

# Import user modules
#from . import Blender_Extract_Skeleton as Skeleton
from cc_blender_plant_volume import blender_plant_rmnoise as rmnoise
from cc_blender_plant_volume import blender_extract_crosssection as section
from cc_blender_plant_volume import blender_extract_attributefiltering as attribute
from cc_blender_plant_volume import blender_point_clustering as cluster
from cc_blender_plant_volume import blender_utility_functions as utility

def cleanup_env(obj_to_remove:list[str] = ["Cube",], type_to_remove:list[str] = ["MESH",]) -> None:
    """Remove non-relevant object from scene"""
    # Unselect all to avoid deleting previously selected object
    utility.select_all(select=False)

    # Loop through object and delete object in the list
    for obj in bpy.data.objects:
        if obj.type in type_to_remove or obj.name in obj_to_remove:
            obj.select_set(True)
            bpy.ops.object.delete(use_global=False)

def import_file(filepath:str, file_ext:str="obj"):
    """Import obj or ply file"""
    # Check if file exist and is obj
    if not os.path.isfile(filepath):
        raise OSError(f"File {filepath} not found")
    if not filepath.endswith(file_ext):
        raise OSError(f"File {filepath} is not an obj")

    # Import file (based on the file_ext)
    import_command = {"ply":"ply_import", "obj":"obj_import"}
    getattr(bpy.ops.wm, import_command[file_ext])(filepath=filepath,
                                                  forward_axis='NEGATIVE_Z',
                                                  up_axis='Y')

def get_location(obj:(bpy.types.Object|None)=None, center_ref:str="surface") -> mathutils.Vector:
    """Reset center and get location of active object or sepified object"""
    # Set the origine reference based on input
    origin_def = {"surface":"ORIGIN_CENTER_OF_MASS",
                  "volume":"ORIGIN_CENTER_OF_VOLUME",}
    # If no input object, work on active object, otherwise, select input object
    if obj is None:
        obj = bpy.context.active_object
    else:
        utility.select_all(select=False)
        bpy.context.view_layer.objects.active = obj
        obj.select_set(True)
    # Set object origin based on pivot point
    bpy.ops.object.origin_set(type=origin_def[center_ref], center='MEDIAN')
    return obj.location

def translate_to_center(ref_obj:bpy.types.Object) -> None:
    """Move all mesh objects to center reference object"""
    bpy.ops.object.mode_set(mode='OBJECT', toggle=False)
    # Set reference object as active
    utility.select_all(select=False)
    bpy.context.view_layer.objects.active = ref_obj
    ref_obj.select_set(True)
    # Get reference object location
    ref_location = get_location()
    # Translate all mesh object based on ref location
    for obj in bpy.data.objects:
        if obj.type == 'MESH':
            obj.location -= ref_location

def get_rotation_to_z(ref_obj:bpy.types.Object, plane_normal:bool=False) -> mathutils.Matrix:
    """Compute rotation matrix vector to allign center of reference object to Z axis"""
    center = get_location(ref_obj)
    if plane_normal:
        # Fit plane to object and allign z axis to plane normal
        plane = cluster.ransac_plane(ref_obj)
        center = mathutils.Vector(plane[0:3])
    # Compute rotation vector (using Euler rotation), then convert to Matrix
    rotation_x =  mathutils.Euler((-np.pi + np.arctan2(center.y, center.z), 0, 0), "XYZ")
    matrix_x   = rotation_x.to_matrix()
    rotation_y = mathutils.Euler((0, np.pi - np.arctan2(center.x, center.z), 0), "XYZ")
    matrix_y   = rotation_y.to_matrix()
    # Transformation by first rotation in Y, then in X
    # (with numpy library @ act as matrix mulitplication operator)
    return matrix_y @ matrix_x

def skew_sym_cross_product(vec1:np.ndarray, vec2:np.ndarray) -> np.ndarray:
    """Return the skew symetric cross product matrix of the cross product of vec1 and vec2"""
    cross = np.cross(vec1, vec2)
    return np.array([[0, -cross[2], cross[1]],
                     [cross[2], 0, -cross[0]],
                     [-cross[1], cross[0], 0]])

# Using https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
def get_rotation_to_axis(ref_obj: bpy.types.Object,
                         target:np.ndarray=np.array([0., 0., 1.])) -> mathutils.Matrix:
    """Compute rotation matrix from normal of reference object to terget vector"""
    # Fit plane to object and allign z axis to plane normal
    ref_vector  = cluster.ransac_plane(ref_obj)[0:3]
    # Compute skew symetric matrix and cos of angle
    skew = skew_sym_cross_product(ref_vector, target)
    cos  = np.dot(ref_vector, target)
    # Return rotation matrix
    return mathutils.Matrix(np.identity(3) + skew + skew**2 / (1+cos))

def apply_rotation(obj:bpy.types.Object, rotation_matrix:mathutils.Matrix) -> None:
    """Apply input Euler rotation vector to object"""
    # Position object around origin for the rotation
    init_location = obj.location.copy()
    obj.location = mathutils.Vector((0, 0, 0))
    obj.rotation_euler.rotate(rotation_matrix)
    # To revert object location, need to multiply the initial location with the rotation matrix
    # (with numpy library @ act as matrix mulitplication operator)
    obj.location = rotation_matrix @ init_location

def allign_to_z(ref_obj:bpy.types.Object, plane_normal=True) ->  int:
    """Rotate object around 2 axis to allign center of ref object to Z axis
    Return -1 if an error happen (otherwise, return 1)
    """
    bpy.ops.object.mode_set(mode='OBJECT', toggle=False)
    # If reference object does not contain any vertices, return -1 and alert the user
    if len(ref_obj.data.vertices) == 0:
        print(f"Warning {ref_obj.name} is empty, cannot be used to for allignment")
        return -1
    # Mark reference object as active
    utility.select_all(select=False)
    bpy.context.view_layer.objects.active = ref_obj
    # Compute rotation matrix differently if using the plane normal or not
    if plane_normal:
        # Apply all previous transformation
        utility.select_all(select=True)
        bpy.ops.object.transform_apply(location=True, rotation=True, scale=True)
        rotation_matrix = get_rotation_to_axis(ref_obj)
        utility.select_all(select=False)
    else:
        # Set pivot point to center of grid
        bpy.context.scene.cursor.location = mathutils.Vector((0, 0, 0))
        bpy.context.scene.tool_settings.transform_pivot_point = 'CURSOR'
        rotation_matrix = get_rotation_to_z(ref_obj)
    # Perform the rotation
    for obj in bpy.data.objects:
        if obj.type == 'MESH':
            apply_rotation(obj, rotation_matrix)
    return 1

def calc_scale_ratio(measure:np.ndarray, expected:float, method:str) -> float:
    """Return scale ratio to scale to the expected value"""
    avg_measure = getattr(np, method, "average")(measure)
    return expected / avg_measure

def scale_to_ref(ratio) -> None:
    """Scale all mesh object to defined ratio"""
    # Initialise transformation
    bpy.ops.object.mode_set(mode='OBJECT', toggle=False)
    utility.select_all(select=False)
    bpy.context.scene.tool_settings.transform_pivot_point = 'CURSOR'
    # Loop through all object and scale mesh objects
    for obj in bpy.data.objects:
        if obj.type == 'MESH':
            # Position object around origin for the rotation
            init_location = obj.location.copy()
            obj.location = mathutils.Vector((0, 0, 0))
            obj.scale = mathutils.Vector((ratio, ratio, ratio))
            obj.location = init_location * ratio
#            obj.select_set(True)
#            bpy.ops.transform.resize(value=(ratio, ratio, ratio), orient_type='GLOBAL',
#                                     center_override=(0.0, 0.0, 0.0))
#            obj.select_set(False)

# Warning: if switch scaling to bmesh, would need to update
def calc_volume(obj:bpy.types.Object) -> float:
    """Return volume from bmesh after application of transformation"""
    # Mark object as active and switch to Edit mode
    bpy.ops.object.mode_set(mode='OBJECT', toggle=False)
    utility.select_all(select=False)
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

def import_model(obj_path:str, file_ext:str="obj") -> None:
    """Prepare the working environment and load the 3D model"""
    # Prepare working environment
    cleanup_env()
    import_file(obj_path, file_ext)

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
    # If no active object alert the user
    assert plant is not None, "No active object"

    # Clean up obj
#    attribute.delete_low_confidence()
    rmnoise.main()
    # Pot and cup detection by skeleton or by color
#    Skeleton.main(name="Pot")
#    Skeleton.delete_isolated()
#    Skeleton.keep_biggest_cluster()
    (pot, cup) = attribute.copy_pot()

    # Allign Object based on pot and cup center
    translate_to_center(ref_obj = pot)
    # Object rotation needs to be aplpied twice, otherwise, allignment not correct (why? to check)
    # Each time call the alignment, record the allignment status
    alignment_status = []
    alignment_status.append(allign_to_z(ref_obj=cup))
    alignment_status.append(allign_to_z(ref_obj=cup))

    # If any of the alignment status returned -1, save error indicator as -1,
    # otherwise, initialise it as 1
    error_indicator = -1 if -1 in alignment_status else 1

    # Compute Cross-section
    cross_section = section.main(plant)
    section_dim = section.get_section_dimension(cross_section)

    # Perform scaling based on measured ratio if the cross-section returned proper values
    if len(section_dim) != 1:
        ratio = calc_scale_ratio(section_dim, pot_size, method="min")
        scale_to_ref(ratio)
    else:
        # section returned -1, skip the scaling and multiply the volume by -1 to alert the user
        error_indicator = -1

    # Create copy of plant excluding pot
    plant_green = attribute.copy_green_plant(plant)
    # Use clusting to only keep biggest cluster
    deleted_vertices = cluster.dbscan_filter(plant_green)
    # If -1 returned as number of deleted vertices, no cluster detected,
    # set error to -1 to indicate error
    if deleted_vertices == -1:
        error_indicator = -1

    # Compute volume
    plant_volume = calc_volume(plant_green) * error_indicator
    print(f"{plant_volume = }")

    # Cleanup unused data
    bpy.ops.outliner.orphans_purge()

    # Return the volume
    return plant_volume

def single_model_prep(obj_path:str,
                      output_path:str="",
                      pot_size:float=0.13,
                      file_ext:str="obj") -> float:
    """Single model preparation: import the model, compute the volume and save output blend file"""
    # Import the model
    import_model(obj_path, file_ext)
    # Prepare the model and compute the volume
    volume = model_prep(pot_size)
    # Save blend file in output folder
    if output_path != "":
        save_blend(obj_path, output_path)
    # Return the volume
    return volume

def process_model_folder(file_list:list[str],
                         working_path:str,
                         volume_path:str,
                         output_path:str="",
                         file_ext:str="obj") -> None:
    """Loop through file list and process all model files"""
    for name in file_list:
        if not name.lower().endswith(file_ext) or "ptscloud" in name.lower():
            continue
        # Process plant model
        model_path = os.path.join(working_path, name)
        volume = single_model_prep(model_path, output_path=output_path, file_ext=file_ext)
        # Save volume to csv file
        with open(volume_path, "a", encoding="utf-8") as volume_file:
            volume_file.write(f"{working_path},{name},{volume}\n")

def loop_through_folders(working_path:str,
                         model_folder:str,
                         model_file_ext:str,
                         output_path:str,
                         output_table:str) -> None:
    """Loop through the folders and process folders matching input name"""
    # Loop through all files and process plant model
    volume_path = os.path.join(output_path, output_table)
    with open(volume_path, "w", encoding="utf-8") as volume_file:
        volume_file.write("Path,Plant name,Volume\n")
    for root, _, files in os.walk(working_path):
        # Check if last folder in root matches with specified model folder
        last_folder = root.split(os.sep)[-1]
        if model_folder == last_folder:
            print(f"Processing {root}")
            process_model_folder(files, root, volume_path, output_path, model_file_ext)

if __name__ == "__main__":
    # Test model preparation on active file
    obj = bpy.context.active_object
    assert obj is not None, "No active object"
    volume = model_prep()
    print(f"{volume=}")
