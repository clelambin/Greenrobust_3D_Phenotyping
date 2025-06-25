"""Model rendering in blender.
As a stand alone, run animation of all model in the scene containing a substring in their names.
"""

# Import libraries
import os
import bpy

# To improve:
# - Add type anotation
# - Create function to loop through object and call visibility or rendering function

# User variables
working_dir = "C:\\Users\\cleme\\Documents\\Hohenheim\\00_Courses\\320_Landscape_and_Plant_Ecology\\MSc_3D_Plant_Characterisation\\3D_Digitalisation\\Rendering"

# User function
def hide_object(obj, hide_bool):
    """Set viewport and render visibility (if boolean set to True, hide, otherwide, show)"""
    obj.hide_set(hide_bool)
    obj.hide_render = hide_bool
    # hide_viewport hide in all viewport, not only in the current one
#    obj.hide_viewport = False

def plant_rendering(obj):
    # Show object
    hide_object(obj, False)
    # Set animation output
    bpy.context.scene.render.filepath = os.path.join(working_dir, obj.name)
    # Start animation
    bpy.ops.render.render(animation=True, use_viewport=True)
    # Hide object
    hide_object(obj, True)

def loop_through_plants(contain_string, call_function, *args, **kwargs):
    """Loop through object containing the defined string and apply the function on each object"""
    # Cannot loop directly through all objects, but itterating through the object list by index worked
    obj_list = bpy.data.collections[0].all_objects
    for obj_index in range(0, len(obj_list)):
        try:
            obj = obj_list[obj_index]
            if contain_string in obj.name:
                print(f"Working on {obj.name}")
                call_function(obj, *args, **kwargs)
            else:
                print(f"Object {obj.name} not including Metashape")
        except IndexError:
            print("No object detected")

if __name__ == "__main__":
    # Stating script
    print("Starting script")
    print("- Initial list:")
    for obj in bpy.data.collections[0].all_objects:
        print(obj.name)
    # hide_view_clear is supposed to reset the view
    #bpy.ops.object.hide_view_clear()

    # Hide all object
    print("- Hide all objects")
    loop_through_plants(contain_string="Metashape", call_function=hide_object, hide_bool=True)

    # Loop through objects, for each, show the object, and start the render
    print("- Start animation rendering")
    loop_through_plants(contain_string="Metashape", call_function=plant_rendering)
