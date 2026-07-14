## Batch thresholding functions

## Libraries
import cv2 as cv                                # Image processing
import os                                       # File manipulation
import numpy as np                              # Array manipulation
from matplotlib import pyplot as plt            # To display multiple images in one plot
import exif                                     # For EXIF manipulation

## Utility functions
def create_folder(folder_name):
    """Check if folder exist, if not, create it, return True if a folder was created"""
    if not os.path.isdir(folder_name):
        os.mkdir(folder_name)
        return True
    else:
        return False

def encode_image(img_array):
    """Encode image array as jpg and returned encoded byte object"""
    # Encode mask as jpg
    _, img_codded = cv.imencode(".jpg", img_array)
    img_bytes     = img_codded.tobytes()
    return(img_bytes)

def find_biggest_contour(contour_list):
    """Compare contour area and output biggest contour"""
    # Extract contour with biggest contour
    biggest_contour = contour_list[0]
    biggest_area    = 0
    for contour in contour_list:
        contour_area = cv.contourArea(contour)
        if  contour_area > biggest_area:
            biggest_contour = contour
            biggest_area    = contour_area
    return(biggest_contour)

## Bath processing class
class batch_threshold:
    """Batch thresholding of images in working folder and generate mask for each processed images"""
    # Setting variables
    img_format = ("jpg")
    input_folder = "jpg"
    # Rotation dictionary: link between EXIF orientation and opencv rotate
    rotate_exif2cv = {3:cv.ROTATE_180,
                      6:cv.ROTATE_90_COUNTERCLOCKWISE,
                      8:cv.ROTATE_90_CLOCKWISE}
    
    def __init__(self, working_folder, mask_folder="jpg_mask", debug=False):
        self.working_folder = working_folder
        self.mask_folder    = mask_folder
        self.input_path = os.path.join(working_folder, self.input_folder)
        self.mask_path  = os.path.join(working_folder, self.mask_folder)

        self.debug=debug
        if debug:
            print(f"Current folder: {os.getcwd()}")
            print(f"Working folder: {self.working_folder}")
            print(f"Input path: {self.input_path}")
            print(f"Mask path: {self.mask_path}")
    
    def init_input_folder(self):
        """Create the input folder and move all images from working folder to input folder"""
        init_working_files = os.listdir(self.working_folder)
        folder_created = create_folder(self.input_path)
        # If folder already exist, input already processed, skipß
        if not folder_created:
            print(f"Folder {self.input_path} already exist, skipping initialisation")
            return
        # Move images from working folder
        for file_name in init_working_files:
            if file_name.lower().endswith(self.img_format):
                os.rename(os.path.join(self.working_folder, file_name), os.path.join(self.input_path, file_name))

    def create_mask(self, img_path, close_dist=5, errode_dist=20):
        """Create mask for current image"""
        ## Load image
        img_bgr  = cv.imread(img_path)
        img_grey = cv.cvtColor(img_bgr, cv.COLOR_BGR2GRAY)

        ## Image Thresholding
        # Apply Otsu threshold
        _, otsu = cv.threshold(img_grey,0,255,cv.THRESH_BINARY+cv.THRESH_OTSU)
        # Apply adaptive threshold
        adapt   = cv.adaptiveThreshold(img_grey,255,cv.ADAPTIVE_THRESH_MEAN_C,cv.THRESH_BINARY,21,2)
        # Combine both mask
        mask = cv.bitwise_and(otsu, adapt)

        ## Mask post-processing
        # Remove noise (apply closing)
        kernel = np.ones((close_dist,close_dist),np.uint8)
        mask = cv.morphologyEx(mask, cv.MORPH_CLOSE, kernel)
        # Apply errosion on mask to increase covered area
        kernel = np.ones((errode_dist,errode_dist),np.uint8)
        mask = cv.erode(mask,kernel,iterations = 1)
        # Invert mask
        mask = cv.bitwise_not(mask)

        # Return generated mask
        return mask

    def rotate_mask(self, img_path, mask):
        """Rotate image based on EXIF orientation value of original image"""
        with open(img_path, "rb") as img_file:
            exif_source = exif.Image(img_file)
        if self.debug and exif_source.has_exif: print("EXIF data found for working image")
        exif_orient = exif_source.get("orientation")
        cv_orient   = self.rotate_exif2cv.get(exif_orient, None)
        if cv_orient is not None:
            if self.debug: print("Rotating mask")
            mask = cv.rotate(mask, cv_orient)
        return mask

    def clip_mask_to_bg(self, mask, close_dist=0, ignore_top_intesect=True):
        """Use findContours to get the white background and clip outer element from initial mask"""
        # Pre-process mask to help detecting white background as contour
        mask_inv    = cv.bitwise_not(mask)
        kernel      = np.ones((close_dist,close_dist),np.uint8)
        mask_close  = cv.morphologyEx(mask_inv, cv.MORPH_CLOSE, kernel)

        # If ignore_top_intesect is true, set the top part of the mask to white to exclude plant from the background
        if ignore_top_intesect:
            (height, width) = mask_close.shape
            # set 20% of the width (centered) and 5% of the height (top) to white
            ignore_x_start = 0
            ignore_x_end   = int(0.5*height)
            ignore_y_start = int((width - 0.2*width)/2)
            ignore_y_end   = int((width + 0.2*width)/2)
            mask_close[ignore_x_start:ignore_x_end, ignore_y_start:ignore_y_end] = 255
        
        # Find biggest contour from pre-processed mask
        contour_list, _ = cv.findContours(mask_close, cv.RETR_EXTERNAL, cv.CHAIN_APPROX_SIMPLE)
        biggest_contour = find_biggest_contour(contour_list)

        # Generate bounding box from contour
        bbox = cv.boundingRect(biggest_contour)
        if self.debug: print(f"Bounding box of white background: {bbox}")

        # Set out of bounding box as hidden (mask = 0)
        (bbox_x, bbox_y, bbox_width, bbox_height) = bbox
        mask[:, 0:bbox_x]    = 0
        mask[:, bbox_x+bbox_width:-1] = 0
        #mask[0:bbox_y, :]    = 0
        #mask[bbox_y+bbox_height:-1, :] = 0

        # Return updated mask
        return mask
        
    # UNUSED: adding orientation to mask is not enough for metashape
    def save_with_exif(self, img_name, img_array, copied_tag=("orientation",)):
        """Copy given list of tag from input image to target image and save new image"""
        with open(os.path.join(self.input_path, img_name), "rb") as img_file:
            exif_source = exif.Image(img_file)
        if self.debug and exif_source.has_exif:
            print("EXIF data found for working image")
        # Encode image array
        img_encoded = encode_image(img_array)
        # Add exif to encoded image
        exif_target = exif.Image(img_encoded)
        # For each tag, saved value to encoded image
        for tag in copied_tag:
            tag_value = exif_source.get(tag)
            exif_target.set(tag, tag_value)
        # Save image in mask folder
        with open(os.path.join(self.mask_path, img_name), "wb") as img_file:
            img_file.write(exif_target.get_file())
    
    def batch_masking(self, clip_bg = False):
        """Generate mask for each image in input folder"""
        folder_created = create_folder(self.mask_path)
        # If folder already exist, alert the user and end the process
        if not folder_created:
            print(f"Folder {self.mask_path} already exist, skipping masking")
            return
        # Process all images in input folder
        for file_name in os.listdir(self.input_path):
            # Skip non-image file
            if not file_name.lower().endswith(self.img_format):
                continue
            # Process image
            img_path = os.path.join(self.input_path, file_name)
            if self.debug: print(f"Working on {img_path}")
            img_mask = self.create_mask(img_path)
            if clip_bg: img_mask = self.clip_mask_to_bg(img_mask)
            img_mask = self.rotate_mask(img_path, img_mask)
            cv.imwrite(os.path.join(self.mask_path, file_name), img_mask)