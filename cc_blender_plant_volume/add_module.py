"""Add module specified in module path to python path
If module not installed in Blender, add_module should be called before calling __init__
"""

import sys

MODULE_PATH = r"/home/clement/Documents/Hohenheim/00_Courses/MSc_3D_Plant_Characterisation/3D_Digitalisation/00_Scripts/MSc_3D_Digitisation_Scripts"
LOCAL_LIBRARY = r"/home/clement/Documents/Hohenheim/00_Courses/MSc_3D_Plant_Characterisation/3D_Digitalisation/00_Scripts/MSc_3D_Digitisation_Scripts/cc_blender_plant_volume/.venv/lib/python3.11/site-packages"


sys.path.append(MODULE_PATH)
sys.path.append(LOCAL_LIBRARY)

def register():
    """Empty register, used for Blender referencing"""

def unregister():
    """Empty unregister, used for Blender referencing"""
