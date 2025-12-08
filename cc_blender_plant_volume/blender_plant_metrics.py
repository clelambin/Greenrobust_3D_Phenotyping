"""Compute metrics relevant to the plant geometry
Measured metrics: - Volume
                  - Total surface area
                  - Dimension (X, Y and Z length of bounding box)
"""

# Import libraries
from typing import Literal
from PIL import Image # Image operations
import bpy            # Blender python
import bmesh          # Blender mesh module
import numpy as np    # Array and matrix operations

# Import user modules
from cc_blender_plant_volume import blender_utility_functions as utility

# Complex type anotation
VIEWPOINT = Literal['LEFT', 'RIGHT', 'BOTTOM', 'TOP', 'FRONT', 'BACK']

def calc_area(mesh:bmesh.types.BMesh) -> float:
    """Return sum of all faces area"""
    return sum(face.calc_area() for face in mesh.faces)

def calc_dimension(mesh: bmesh.types.BMesh, ref:np.ndarray|None = None) -> np.ndarray:
    """Return X, Y and Z dimension of the bounding box of input mesh"""
    # Compute min and max for all 3 dimenions of the vertices coordinate
    vertex_coord = np.array([vertex.co for vertex in mesh.verts])
    # If no reference point specified, compare to min dimension
    if ref is None:
        ref = np.min(vertex_coord, axis=0)
    return np.max(vertex_coord, axis=0) - ref

# From https://blender.stackexchange.com/questions/290437/set-orthographic-top-view-with-python
# Warning: - Relying heavily on operation (context dependant, less stable, ...)
def uv_project(obj:bpy.types.Object,
               tmp_path:str,
               projection_view:VIEWPOINT="TOP") -> None:
    """Create UV projection of input object from projection viewpoint"""
    # Edit input object
    utility.make_active(obj)
    view3d = utility.get_view3d()
    # Set viewpoint to projected view and project from view within correct context
    # (need to use with to tell blender from which area (the view3d) to work)
    # (the last regions is the 3d windows)
    with bpy.context.temp_override(area=view3d, region=view3d.regions[-1]):
        bpy.ops.view3d.view_axis(type=projection_view)
    # Need to update the view3d, otherwise, viewpoint not considered for uv_project_from_view
    view3d.spaces[0].region_3d.update()
    bpy.ops.object.mode_set(mode='EDIT', toggle=False)
    with bpy.context.temp_override(area=view3d, region=view3d.regions[-1]):
        bpy.ops.mesh.select_all(action="SELECT")
        bpy.ops.uv.project_from_view(orthographic=True,
                                     camera_bounds=False,
                                     correct_aspect=True,
                                     scale_to_bounds=True)
    # Export UV map as image
    bpy.ops.uv.export_layout(filepath=tmp_path, size=(480, 480), opacity=1)
    # Switch back to objetc mode
    bpy.ops.object.mode_set(mode='OBJECT', toggle=False)

# (image type in PIL library confusing: PIL.Image = Module, PIL.Image.Image = image type)
def image_pixel_ratio(img:Image.Image) -> float:
    """Return the ratio of non transparent pixel in input image"""
    # Extract only alpha channel from image
    img_alpha = img.getchannel("A")
    # Load alpha image and convert to numpy
    img_alpha.load()
    img_array = np.asarray(img_alpha)
    img_alpha.close()
    # Return ration of non-zero alpha pixels
    return np.count_nonzero(img_array) / img_array.size

def area_project(obj: bpy.types.Object,
                 tmp_path:str,
                 projection_view:VIEWPOINT="TOP") -> float:
    """Use UV to project object onto a projection view and cound non-transparent pixel ratio"""
    # Create UV projection of object into given viewpoint
    uv_project(obj, tmp_path, projection_view)
    # Open image and return non-transparent ratio
    img = Image.open(tmp_path)
    return image_pixel_ratio(img)

# Warning: if switch scaling to bmesh, would need to update
def calc_metrics(obj:bpy.types.Object|None=None,
                 img_path:str|None=None,
                 ref:np.ndarray|None=None) -> dict[str, float]:
    """Return metrics (volume, area, dimension) from bmesh after application of transformation
    Return metrics are saved within a dictionary
    """
    # Initialise metrics dictionary
    metrics = {
        "Volume"         : 0.,
        "Dim_X"          : 0.,
        "Dim_Y"          : 0.,
        "Dim_Z"          : 0.,
        "Height"         : 0.,
        "Cumul_Area"     : 0.,
        "Font_Area_Ratio": 0.,
        "Left_Area_Ratio": 0.,
        "Top_Area_Ratio" : 0.,
        "Code_Error"     : 0,
    }
    # If no object, return empty dictionary (for the keys)
    if obj is None:
        return metrics
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
    # Calculate mesh metrics
    metrics["Volume"] = temp_mesh.calc_volume()
    metrics["Cumul_Area"] = calc_area(temp_mesh)
    metrics["Dim_X"], metrics["Dim_Y"], metrics["Dim_Z"] = calc_dimension(temp_mesh)
    _, _, metrics["Height"] = calc_dimension(temp_mesh, ref=ref)
    # Free up bmesh and switch back to Object mode
    orig_mesh.free()
    temp_mesh.free()
    bpy.ops.object.mode_set(mode='OBJECT', toggle=False)
    # Calculate projection metrics
    project_view:list[VIEWPOINT] = ["FRONT", "LEFT", "TOP"]
    project_area:list[float] = [0.] * len(project_view)
    if img_path is not None:
        for index, view in enumerate(project_view):
            project_area[index] = area_project(obj, img_path, view)
    metrics["Font_Area_Ratio"], metrics["Left_Area_Ratio"], metrics["Top_Area_Ratio"] = project_area
    # Return all computed metrix
    return metrics

if __name__ == "__main__":
    IMG_PATH = r"C:\Users\cleme\Documents\Hohenheim\00_Courses\320_Landscape_and_Plant_Ecology\MSc_3D_Plant_Characterisation\3D_Digitalisation\00_Scripts\TMP_Blender_Allignment_Issue\test.png"
    test_obj = bpy.context.active_object
    assert test_obj is not None, "No active object"
    print(f"{calc_metrics(test_obj, IMG_PATH)=}")
#    # Test projected area only
#    print(f"{area_project(test_obj, IMG_PATH)=}")
