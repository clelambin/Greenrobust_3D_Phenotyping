"""Model rendering in blender.
As a stand alone, run animation of all model in the scene containing a substring in their names.
"""

# Import libraries
from collections.abc import Callable
from math import pi
import os
import bpy

# Import user modules
from cc_blender_plant_volume.blender_utility_functions import hide_object

# To improve:
# - Add type anotation
# - Create function to loop through object and call visibility or rendering function

# User variables
WORKING_DIR = r"C:\Users\cleme\Documents\Hohenheim\00_Courses\320_Landscape_and_Plant_Ecology\MSc_3D_Plant_Characterisation\3D_Digitalisation\Rendering"

# User function
def plant_rendering(obj:bpy.types.Object,
                    output_dir:str,
                    render_animation:bool=True):
    """Create rendering of input object"""
    # Show object
    hide_object(obj, False)
    # Set animation output
    bpy.context.scene.render.filepath = os.path.join(output_dir, obj.name)
    # Start animation
    # (write_still save the image to output folder if not an animation,
    #  for animation, render always saved)
    bpy.ops.render.render(animation=render_animation,
                          use_viewport=True,
                          write_still=True)
    # Hide object
    hide_object(obj, True)

def is_target_object(obj: bpy.types.Object,
                     contain_string: str|None,
                     target_type: str|None) -> bool:
    """Return True if object name contain inut string or is of input type
    If both criteria defined, both used in criteria
    """
    # Both criteria cannot be None at the same time
    assert contain_string is not None or target_type is not None, "No target criteria"
    # If either criteria not defined, used the other one as criteria
    if contain_string is None:
        return obj.type == target_type
    if target_type is None:
        return contain_string in obj.name
    # Otherwise, use both criteria to check if object should be processed
    return obj.type == target_type and contain_string in obj.name

def loop_through_plants(call_function:Callable,
                        *args,
                        contain_string:str|None=None,
                        target_type:str|None=None,
                        **kwargs):
    """Loop through object containing the defined string and apply the function on each object"""
    for obj in bpy.data.objects:
        if is_target_object(obj, contain_string, target_type):
            print(f"Working on {obj.name}")
            call_function(obj, *args, **kwargs)
        else:
            print(f"Object {obj.name} not including {contain_string}")

def material_from_attribute(material_name:str="ColorAttribute", attribute_name:str="Col") -> bpy.types.Material:
    """Create material from attribute using Shadder Node"""
    # Initialise material node
    material = bpy.data.materials.new(name = material_name)
    material.use_nodes = True
    material_node = material.node_tree.nodes
    # Node Principled BSDF
    principled_bsdf = material_node.new("ShaderNodeBsdfPrincipled")
    principled_bsdf.name = "Principled BSDF"
    principled_bsdf.distribution = 'MULTI_GGX'
    principled_bsdf.subsurface_method = 'RANDOM_WALK'
    # Node Material Output
    material_output = material_node.new("ShaderNodeOutputMaterial")
    material_output.name = "Material Output"
    material_output.is_active_output = True
    material_output.target = 'ALL'
    # Node Attribute
    attribute = material_node.new("ShaderNodeAttribute")
    attribute.name = "Attribute"
    attribute.attribute_name = attribute_name
    attribute.attribute_type = 'GEOMETRY'
    # Initialize material links
    material.node_tree.links.new(principled_bsdf.outputs[0], material_output.inputs[0])
    material.node_tree.links.new(attribute.outputs[0], principled_bsdf.inputs[0])
    # Return created material
    return material

def assign_material(obj:bpy.types.Object, material:bpy.types.Material) -> None:
    """Assign material to input object"""
    material_slot = obj.data.materials.find(material.name)
    # If material not yet assigned, assign it to input object
    if material_slot == -1:
        obj.data.materials.append(material)
    # Update material slot to match the last index
        material_slot = len(obj.data.materials) - 1
    # Make assigned material as active
    obj.active_material_index = material_slot

def set_material_from_attribute(material_name:str="ColorAttribute") -> None:
    """Loop through existing mesh and assign the material based on color attribute"""
    # Create the material if it does not exists
    material_list = bpy.data.materials
    if material_list.find(material_name) == -1:
        material = material_from_attribute(material_name)
    else:
        material = material_list[material_name]
    # Assign material to all meshs in the scene
    loop_through_plants(assign_material, target_type="MESH", material=material)

def model_rendering(output_dir:str, scene_name:str) -> None:
    """Create rendering of current prepared model"""
    # Initialise scene object
    scene = bpy.context.scene
    assert scene is not None, "No active scene found"
    camera = scene.camera
    assert camera is not None, "No Camera in the scene"
    # Set background color to white and image resolution
    bpy.data.worlds["World"].node_tree.nodes["Background"].inputs[0].default_value = (1, 1, 1, 1)
    bpy.context.scene.render.resolution_x = 360
    bpy.context.scene.render.resolution_y = 360
    # Hide initial mesh from the scene
    loop_through_plants(call_function=hide_object,
                        contain_string="Metashape",
                        target_type="MESH",
                        hide_bool=True)
    # Assign mmaterial from color attribute to all Mesh
    set_material_from_attribute()
    # Define camera viewpoints
    camera_vp = {
         "front":{"location":(0, -2, 0.5),
                  "rotation":(pi/2, 0, 0)},
         "left" :{"location":(-2, 0, 0.5),
                  "rotation":(pi/2, 0, -pi/2)},
         "top"  :{"location":(0, 0, 1.5),
                  "rotation":(0, 0, 0)},
    }
    # Create animation for each camera viewpoints
    for view, position in camera_vp.items():
        # Define output path
        img_name = f"{scene_name}_{view}.jpg"
        scene.render.filepath = os.path.join(output_dir, img_name)
        # Set camera position and save rendering
        camera.location = position["location"]
        camera.rotation_euler = position["rotation"]
        bpy.ops.render.render(write_still=True, use_viewport=True)
    # Reset model visibility
    loop_through_plants(call_function=hide_object,
                        contain_string="Metashape",
                        target_type="MESH",
                        hide_bool=False)

def main():
    """Loop through object and start rendering making only current object visible"""
    # Stating script
    print("Starting script")
    print("- Initial list:")
    for obj in bpy.data.objects:
        print(obj.name)
    # hide_view_clear is supposed to reset the view
    #bpy.ops.object.hide_view_clear()

    # Hide all object
    print("- Hide all objects")
    loop_through_plants(call_function=hide_object,
                        contain_string="Metashape",
                        hide_bool=True)

    # Loop through objects, for each, show the object, and start the render
    print("- Start animation rendering")
    loop_through_plants(call_function=plant_rendering,
                        contain_string="Metashape",
                        output_dir=WORKING_DIR)

if __name__ == "__main__":
    main()
