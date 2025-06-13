## Metashape batchprocessing
# Metashape script to prepare 3D digitisation from different image set

# User variables
#working_dir = r"D:\Clement\3D_Digitalisation\2025-05-28_Harvest_Solanum_dulcamara"
working_dir = r"."
output_dir  = r"02_Metashape_BatchProc_withConfidence"

# Libraries
import Metashape
import os
import re
from datetime import datetime

# Script function
def create_folder(path):
    """Create given folder if it does not already exists"""
    if not os.path.isdir(path):
        os.makedirs(path)

## Metashape processing class
class metashape_proc:

    # Setting variables
    img_format = ("jpg")
    # List of possible folder to process: the last folder in the list has precedence
    img_possible_folder  = ("jpg", "jpg_filtered")
    mask_possible_folder = ("jpg_mask", "jpg_mask_erode")
    
    def __init__(self, working_dir=".", output_dir=".", output_tag="Metashape_Output", new_env=True):
        """Initialise the Metashape document"""

        # Access current document and clear it (equivalent to File>New)
        self.doc = Metashape.app.document
        if new_env:
            self.doc.clear()
            # Create new chunk
            self.chunk = self.doc.addChunk()
        else:
            # Work on existing chunk
            self.chunk = self.doc.chunk
        # Set directories
        self.active_dir  = os.getcwd()
        self.working_dir = working_dir
        self.output_dir  = output_dir
        # Set output tag
        date_tag = datetime.today().strftime('%Y%m%d_%H%M')
        self.output_tag = f"{output_tag}_{date_tag}"
        # Check which folder to use as input
        self.init_input_folder()

    def init_input_folder(self):
        """Setup input folders used to import photos"""
        # Convert folder name to relative path from active directory
        img_possible_path  = [os.path.join(self.working_dir, folder_name) for folder_name in self.img_possible_folder]
        mask_possible_path = [os.path.join(self.working_dir, folder_name) for folder_name in self.mask_possible_folder]
        # Look for existing folder to process
        self.img_folder = ""
        self.mask_folder = ""
        for folder in img_possible_path:
            if os.path.isdir(folder):
                self.img_folder = folder
        for folder in mask_possible_path:
            if os.path.isdir(folder):
                self.mask_folder = folder
        # If no image list found, raise an error, if no mask found, alert the user
        if self.img_folder == "":
            raise NameError(f"No image folder found matching one of the following: {self.img_possible_folder}")
        if self.mask_folder == "":
            print(f"No mask folder found matching one of the following: {self.mask_possible_folder}")

    def import_photos(self):
        """Import list of images to chunk and corresponding mask"""
        file_list = os.listdir(self.img_folder)
        img_list = [file_name for file_name in file_list if file_name.lower().endswith(self.img_format)]
        img_path = [os.path.join(self.img_folder, img_name) for img_name in img_list]
        self.chunk.addPhotos(img_path)
        if os.path.isdir(self.mask_folder):
            # move to mask folder temporarily to import the mask
            os.chdir(self.mask_folder)
            self.chunk.generateMasks(path="{filename}.JPG", masking_mode=Metashape.MaskingMode.MaskingModeFile)
            os.chdir(self.active_dir)

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

    def model_from_pictures(self):
        """Full model preparation on the loaded pictures"""
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
            # Save project both as archive (psz) for transfering to another machine and open (psx) for editing
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
    for file_name in file_list:
        if os.path.isdir(file_name) and file_name.startswith(("BR", "HS", "HV", "NB", "SD", "SL", "TA")):
            print(f"Working on {file_name}")
            # Extract plant code from file name
            try:
                plant_code = re.match("[A-Z]{2}[0-9]{3}_[A-Z][0-9]{2}", file_name)[0]
            except TypeError:
                # No match found, use full folder name instead
                plant_code = file_name
            current_process = metashape_proc(working_dir=file_name,
                                             output_dir=output_dir,
                                             output_tag=f"Metashape_{plant_code}")
            current_process.import_photos()
            current_process.model_from_pictures()
            current_process.save_output()