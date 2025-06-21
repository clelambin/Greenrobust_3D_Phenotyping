"""Add module specified in module path to python path
If module not installed in Blender, add_module should be called before calling __init__
"""

import sys

MODULE_PATH = r"C:\Users\cleme\Documents\Hohenheim\00_Courses\320_Landscape_and_Plant_Ecology\MSc_3D_Plant_Characterisation\3D_Digitalisation\00_Scripts\Mcs_3D_Digitisation_Scripts"
LOCAL_LIBRARY = r"C:\Users\cleme\AppData\Roaming\Python\Python311\site-packages"

sys.path.append(MODULE_PATH)
sys.path.append(LOCAL_LIBRARY)

def register():
    pass

def unregister():
    pass
