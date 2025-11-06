## Metashape batchprocessing
# Metashape script to prepare 3D digitisation from different image set

# User variables
#working_dir = r"D:\Clement\3D_Digitalisation\2025-05-28_Harvest_Solanum_dulcamara"
working_dir = r"."
output_dir  = r"03_Metashape_BatchProc_withMarkers"

# Libraries
import Metashape
import os
import re
import math
from datetime import datetime

# Script function
def create_folder(path):
    """Create given folder if it does not already exists"""
    if not os.path.isdir(path):
        os.makedirs(path)

def cart2cyl(x:float, y:float, z:float=0):
    """Convert cartesian coordinate to cylindrical coordinate"""
    rho = math.sqrt(x**2 + y**2)
    phi = math.atan2(y, x)
    return(rho, phi, z)

def cyl2cart(rho:float, phi:float, z:float=0):
    """Convert cylindrical coordinate to cartesian coordinate"""
    x = rho * math.cos(phi)
    y = rho * math.sin(phi)
    return(x, y, z)

## Metashape processing class
class metashape_proc:

    # Setting variables
    img_format = ("jpg")
    # List of possible folder to process: the last folder in the list has precedence
    img_possible_folder  = ("jpg", "jpg_filtered", "Camera_high", "Camera_mid", "Camera_low")
    mask_possible_folder = ("jpg_mask", "jpg_mask_erode")
    
    def __init__(self, working_dir=".", output_dir=".", output_tag="Metashape_BatchProc", new_env=True):
        """Initialise the Metashape document"""

        # Access current document and clear it (equivalent to File>New)
        self.doc = Metashape.app.document
        if new_env:
            self.doc.clear()
#            # Create new chunk
#            self.chunk = self.doc.addChunk()
#        else:
#            # Work on existing chunk
#            self.chunk = self.doc.chunk
        # Initialise import chunk
        self.import_chunk = []
        # Set directories
        self.active_dir  = os.getcwd()
        self.working_dir = working_dir
        self.output_dir  = output_dir
        # Set output tag
        date_tag = datetime.today().strftime('%Y%m%d_%H%M')
        self.output_tag = f"{output_tag}_{date_tag}"
        # Save log into working folder
        Metashape.app.settings.log_enable = True
        Metashape.app.settings.log_path = os.path.join(output_dir, f"{output_tag}_Log_{date_tag}.txt")
        # Check which folder to use as input
        self.init_input_folder()

    def init_input_folder(self):
        """Setup input folders used to import photos"""
        # Convert folder name to relative path from active directory
        img_possible_path  = [os.path.join(self.working_dir, folder_name) for folder_name in self.img_possible_folder]
        mask_possible_path = [os.path.join(self.working_dir, folder_name) for folder_name in self.mask_possible_folder]
        # Look for existing folder to process
        self.img_folders = []
        self.mask_folders = []
        for folder in img_possible_path:
            if os.path.isdir(folder):
                self.img_folders.append(folder)
        for folder in mask_possible_path:
            if os.path.isdir(folder):
                self.mask_folders.append(folder)
        # If no image list found, raise an error, if no mask found, attempt to generate a mask based on image background
        if len(self.img_folders) == 0:
            raise NameError(f"No image folder found matching one of the following: {self.img_possible_folder}")
        if len(self.mask_folders) == 0:
            print(f"No mask folder found matching one of the following: {self.mask_possible_folder}")

    def import_photos(self):
        """Import list of images to chunk and corresponding mask"""
        print(self.img_folders)
        for img_folder in self.img_folders:
            # Create chunk for each import folder
            self.import_chunk.append(self.doc.addChunk())
            file_list = os.listdir(img_folder)
            img_list = [file_name for file_name in file_list if file_name.lower().endswith(self.img_format)]
            img_path = [os.path.join(img_folder, img_name) for img_name in img_list]
            self.import_chunk[-1].addPhotos(img_path)

            # Look for corresponding background image and generate diff mask
            img_background = os.path.basename(f"{img_folder}_background.JPG")
            if os.path.isfile(img_background):
                self.import_chunk[-1].generateMasks(path=img_background,
                                                    masking_mode=Metashape.MaskingMode.MaskingModeBackground,
                                                    tolerance=10)
            else:
                print(f"No background found at path: {img_background}")

        # Merge all import chunks and set merged chunk as active
        self.doc.mergeChunks(chunks=self.import_chunk, copy_masks=True, merge_assets=True)
        self.doc.chunk = self.doc.chunks[-1]
        self.chunk = self.doc.chunk

        # Apply masks from mask folder on remaining images without mask
        camera_set = set(self.chunk.cameras)
        if self.chunk.masks is None:
            # - if no masks found, set the list of camera without mask to the whole list of camera
            camera_nomask = list(camera_set)
        else:
            # - look for images without mask in merged chunk
            mask_set = set(self.chunk.masks.keys())
            camera_nomask = list(camera_set.difference(mask_set))
        
        # - apply mask (if any) on given images
        for mask_folder in self.mask_folders:
            # move to mask folder temporarily to import the mask
            os.chdir(mask_folder)
            self.chunk.generateMasks(path="{filename}.JPG",
                                     masking_mode=Metashape.MaskingMode.MaskingModeFile,
                                     cameras=camera_nomask)
            os.chdir(self.active_dir)

    def marker_ref_circular(self, radius:float, angle:float, start_index:int=0, start_angle:float=0.0) -> int:
        """Update marker reference location to fit in a circle with polar coordinate (radius, angle)"""
        # Initialise marker angle
        marker_angle:float = start_angle
        marker_index:int   = start_index
        while marker_angle-start_angle < 2*math.pi:
            # Compute marker location in cartesian coordinate
            marker_polar = (radius, marker_angle)
            marker_coord = Metashape.Vector(cyl2cart(*marker_polar))
            # Update marker location
            self.chunk.markers[marker_index].reference.location = marker_coord
            self.chunk.markers[marker_index].reference.enabled = True
            # Increment marker angle
            marker_angle += angle
            marker_index += 1
        # Return next (non-updated) marker
        return marker_index

    def detect_markers(self, target_type:str="CircularTarget12bit"):
        """Detect markers and assign coordinates"""
        # Detect marker
        marker_target = getattr(Metashape.TargetType, target_type)
        self.chunk.detectMarkers(target_type=marker_target, tolerance=50)
        # Update marker reference location
        # - The 12 first markers are equi distance in a circule of radius 0.18m separated by an angle of pi/6 rad
        next_marker = self.marker_ref_circular(radius=0.18, angle=math.pi/6)
#        # - The next 4 markers are spaced in a square of diagonal 0.2m
#        self.marker_ref_circular(radius=0.1, angle=math.pi/2, start_index=next_marker, start_angle=math.pi/2)

    def align_photos(self):
        """Allign photo operation"""
        self.chunk.matchPhotos(downscale = 1,            # High image alignment accuracy
                               generic_preselection = True,
                               reference_preselection = False,
                               filter_mask = True,
                               mask_tiepoints = False,
                               filter_stationary_points = True,
                               keypoint_limit = 40000,
                               keypoint_limit_per_mpx = 1000,
                               tiepoint_limit = 4000,
                               guided_matching = False,
                               reset_matches=True)
        self.chunk.alignCameras(adaptive_fitting = True)

    def build_model(self):
        """Build model operation"""
        self.chunk.buildDepthMaps(downscale = 2,         # High depth map quality
                                  filter_mode = Metashape.FilterMode.MildFiltering,
                                  reuse_depth = False)
        self.chunk.buildModel(surface_type = Metashape.SurfaceType.Arbitrary,
                              interpolation = Metashape.Interpolation.EnabledInterpolation,
                              face_count = Metashape.FaceCount.MediumFaceCount,
                              source_data = Metashape.DataSource.DepthMapsData,
                              vertex_colors = True,
                              vertex_confidence = True,
                              volumetric_masks = False)

    def build_texture(self):
        """Build texture operation"""
        self.chunk.buildUV(mapping_mode = Metashape.GenericMapping)
        self.chunk.buildTexture(blending_mode = Metashape.BlendingMode.MosaicBlending,
                               texture_size = 8192,
                               fill_holes = True,
                               ghosting_filter = True,
                               texture_type = Metashape.Model.TextureType.DiffuseMap,
                               transfer_texture = False)

    def model_from_pictures(self, markers=True):
        """Full model preparation on the loaded pictures"""
        if markers: self.detect_markers()
        self.align_photos()
        self.build_model()
        self.build_texture()

    def save_output(self, save_project=True, export_model=True):
        """Export generated model and save metashape project"""
        if export_model:
            # Save model both in ply and obj
            # ply contain cofidence information
            # obj contain texture information
            model_file_ply = f"{self.output_tag}.ply"
            model_file_obj = f"{self.output_tag}.obj"
            self.chunk.exportModel(path=os.path.join(self.output_dir, model_file_ply),
                                   save_confidence=True, save_colors=True)
            self.chunk.exportModel(path=os.path.join(self.output_dir, model_file_obj),
                                   save_confidence=True, save_colors=True)
        if save_project:
            # Save project both as archive (psz) for data transfer and open (psx) for editing
            project_file_psx = f"{self.output_tag}.psx"
            project_file_psz = f"{self.output_tag}.psz"
            self.doc.save(os.path.join(self.output_dir, project_file_psx))
            self.doc.save(os.path.join(self.output_dir, project_file_psz))

# Code block: run if script called directly)
if __name__ == "__main__":
    # Move to yorking directory
    os.chdir(working_dir)

    # Create ouput directory if it does not exists
    create_folder(output_dir)

    # Processing all folder within working directory
    file_list = os.listdir(".")
    output_list = os.listdir(output_dir)
    for file_name in file_list:
        if os.path.isdir(file_name) and file_name.startswith(("AT", "BR", "HS", "HV", "NB", "SD", "SL", "TA")):
            print(f"Working on {file_name}")
            # Extract plant code from file name
            plant_find = re.match("[A-Z]{2}[0-9]{3}_[A-Z][0-9]{2}", file_name)
            # If no match found, use full folder name, otherwise, use first match
            plant_code = file_name if plant_find is None else plant_find[0]
            # Check if plant code already process, if so, skip it
            plant_check = re.compile(f"Metashape_{plant_code}.*")
            if len(list(filter(plant_check.match, output_list))) > 0:
                print(f"Plant {plant_code} already processed, skipping")
                continue
            current_process = metashape_proc(working_dir=file_name,
                                             output_dir=output_dir,
                                             output_tag=f"Metashape_{plant_code}")
            current_process.import_photos()
            current_process.model_from_pictures()
            current_process.save_output()
