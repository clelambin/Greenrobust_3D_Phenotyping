"""Iterative face collapse to derive the rough skeleton from the 3D model.
Initially used to detect pot location.
(Currently replaced by color attribute to detect the different componants)
"""

import bpy
import numpy as np

# To do:
# - Convert to bmesh to speed up process

def face_collapse(select_ratio:float = 0.5, collapse_dist:float = 0.01) -> None:
    """Randomly select faces and apply remove double on them"""
    # Switch to edit mode
    bpy.ops.object.mode_set(mode='EDIT', toggle=False)
    # Unselect all
    bpy.ops.mesh.select_all(action='DESELECT')
    # Randomly select a portion of the faces to apply the remove double
    bpy.ops.mesh.select_mode(type='FACE')
    bpy.ops.mesh.select_random(ratio=select_ratio)
    bpy.ops.mesh.remove_doubles(threshold=collapse_dist)

def hide_loose():
    """Hide loose faces from active view"""
    # Switch to edit mode
    bpy.ops.object.mode_set(mode='EDIT', toggle=False)
    # Unselect all
    bpy.ops.mesh.select_all(action='DESELECT')
    # Select loose faces
    bpy.ops.mesh.select_mode(type='FACE')
    bpy.ops.mesh.select_loose()
    # Hide selected
    bpy.ops.mesh.hide(unselected=False)

def itterative_collapse(min_collapse:float = 0.01, collapse_increment:float = 0.01, nb_itteration:int = 10) -> None:
    """Itteratively run face_collapse, increasing the collapse distance"""
    for itter in range(nb_itteration):
        collapse_dist = min_collapse + collapse_increment*itter
        face_collapse(collapse_dist=collapse_dist)
        hide_loose()

def delete_isolated() -> None:
    """Delete isolated faces and vertices not linked to a face"""
    # Switch to edit mode
    bpy.ops.object.mode_set(mode='EDIT', toggle=False)
    # Hide unconnected faces
    bpy.ops.mesh.select_all(action='DESELECT')
    hide_loose()
    # Select all remaining faces and hide unselected
    bpy.ops.mesh.select_mode(type='FACE')
    bpy.ops.mesh.select_random(ratio=1)
    bpy.ops.mesh.select_mode(type='VERT')
    bpy.ops.mesh.hide(unselected=True)
    # Show and select all hidden, then delete them
    bpy.ops.mesh.select_all(action='DESELECT')
    bpy.ops.mesh.reveal(select=True)
    bpy.ops.mesh.delete(type='VERT')
    # Switch back to object mode
    bpy.ops.object.mode_set(mode='OBJECT', toggle=False)

def keep_biggest_cluster(obj:(bpy.types.Object|None)=None):
    """Select biggest cluster by looking at points closest to median coordinate, then delete remainnig"""
    # If no object given in input, work on active object, otherwise, set input object as active
    if obj is None:
        obj = bpy.context.active_object
    else:
        bpy.ops.object.select_all(action='DESELECT')
        bpy.context.view_layer.objects.active = obj
        obj.select_set(True)
    # Switch to edit mode
    bpy.ops.object.mode_set(mode='EDIT', toggle=False)
    # Unselect all vertices
    bpy.ops.mesh.select_all(action='DESELECT')
    # Extract vertex median coordinates
    vertex_coord  = np.array([vertex.co for vertex in obj.data.vertices])
    vertex_median = np.median(vertex_coord, axis=0)
    # Select closest vertex to median coordinate
    dist_to_med   = np.linalg.norm(vertex_coord - vertex_median, axis=1)
    # Need to be in Object mode to select vertex
    bpy.ops.object.mode_set(mode='OBJECT', toggle=False)
    obj.data.vertices[np.argmin(dist_to_med)].select = True
    bpy.ops.object.mode_set(mode='EDIT', toggle=False)
    bpy.ops.mesh.select_linked()
    # Inverse selection and delete
    bpy.ops.mesh.select_all(action='INVERT')
    bpy.ops.mesh.delete(type='VERT')
    # Switch back to object mode
    bpy.ops.object.mode_set(mode='OBJECT', toggle=False)

def main(name:str = "Skeleton") -> None:
    """Main function, generate the skeleton from active object"""
    # Duplicate object
    bpy.ops.object.mode_set(mode='OBJECT', toggle=False)
    bpy.ops.object.duplicate()

    # Itteratively collapse the mesh to get the skeleton
    itterative_collapse(min_collapse=0.02,
                        collapse_increment=0.02,
                        nb_itteration=20)

    # Rename duplicated object
    bpy.ops.object.mode_set(mode='OBJECT', toggle=False)
    dup_obj = bpy.context.active_object
    dup_obj.name = name

if __name__ == "__main__":
    main(name="Pot")
    delete_isolated()
    keep_biggest_cluster()