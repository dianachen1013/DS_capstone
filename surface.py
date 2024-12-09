from nilearn import plotting, surface, datasets
import nibabel as nib
import numpy as np
import os
from matplotlib.colors import ListedColormap

# Path to Schaefer atlas
schaefer_file = 'D:/ibp/schaefer/Schaefer2018_400Parcels_7Networks_order_FSLMNI152_1mm.nii/Schaefer2018_400Parcels_7Networks_order_FSLMNI152_1mm.nii'

# Load the Schaefer atlas
atlas_img = nib.load(schaefer_file)

# Define clusters
cluster_rois = {
    "Cluster 1": [1, 36, 38, 39, 49, 68, 69, 77, 85, 86, 87, 88, 89, 90, 92, 94, 95, 96, 102, 103, 104, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 126, 127, 131, 132, 134, 135, 137, 139, 140, 141, 147, 148, 149, 150, 153, 154, 157, 158, 159, 161, 162, 165, 166, 167, 169, 171, 174, 177, 178, 180, 181, 182, 184, 187, 195, 232, 238, 242, 271, 274, 290, 295, 296, 297, 299, 300, 303, 304, 305, 308, 309, 319, 321, 322, 323, 326, 327, 328, 330, 331, 332, 333, 337, 338, 339, 340, 342, 343, 344, 346, 347, 349, 350, 353, 354, 359, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 373, 374, 375, 376, 377, 378, 380, 382, 383, 385, 386, 387, 388, 389, 390, 393, 399],
    "Cluster 2": [2, 3, 4, 5, 7, 8, 10, 11, 13, 14, 16, 17, 18, 22, 23, 26, 27, 28, 30, 41, 51, 52, 56, 58, 60, 70, 71, 72, 73, 74, 75, 78, 79, 80, 81, 93, 138, 145, 160, 168, 170, 172, 175, 179, 183, 201, 203, 204, 205, 206, 208, 210, 211, 213, 215, 218, 221, 222, 223, 226, 227, 229, 275, 280, 282, 284],
    "Cluster 3": [33, 34, 35, 37, 40, 42, 43, 44, 45, 46, 47, 48, 50, 53, 54, 55, 57, 59, 61, 62, 63, 64, 65, 66, 67, 83, 84, 91, 100, 101, 105, 111, 136, 142, 151, 152, 155, 156, 188, 209, 216, 228, 230, 235, 236, 237, 239, 240, 241, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255, 256, 257, 258, 259, 260, 261, 262, 263, 264, 265, 266, 267, 268, 269, 270, 272, 273, 276, 278, 279, 283, 286, 287, 288, 289, 293, 298, 301, 306, 307, 310, 312, 314, 315, 316, 317, 318, 335, 341, 348, 351, 352],
    "Cluster 4": [6, 9, 12, 15, 19, 20, 21, 24, 25, 29, 31, 32, 76, 82, 97, 98, 99, 106, 107, 108, 109, 110, 125, 128, 129, 130, 133, 143, 144, 146, 163, 164, 173, 176, 185, 186, 189, 190, 191, 192, 193, 194, 196, 197, 198, 199, 200, 202, 207, 212, 214, 217, 219, 220, 224, 225, 231, 233, 234, 277, 281, 285, 291, 292, 294, 302, 311, 313, 320, 324, 325, 329, 334, 336, 345, 355, 356, 357, 358, 360, 361, 379, 381, 384, 391, 392, 394, 395, 396, 397, 398, 400],
}

# Create a new data array to represent clusters as integers
atlas_data = atlas_img.get_fdata()
highlighted_data = np.zeros_like(atlas_data)

# Assign cluster values
for cluster_idx, rois in enumerate(cluster_rois.values(), start=1):
    for roi in rois:
        highlighted_data[atlas_data == roi] = cluster_idx

# Save the updated atlas
highlighted_img = nib.Nifti1Image(highlighted_data, affine=atlas_img.affine)
highlighted_atlas_path = os.path.join(os.getcwd(), "highlighted_atlas.nii.gz")
nib.save(highlighted_img, highlighted_atlas_path)

# Fetch fsaverage mesh for plotting
fsaverage = datasets.fetch_surf_fsaverage()

# Project volume to surface for both hemispheres
texture_left = surface.vol_to_surf(highlighted_img, fsaverage.pial_left)
texture_right = surface.vol_to_surf(highlighted_img, fsaverage.pial_right)

# Define the custom discrete colormap
cmap = ListedColormap(["#E91E63", "#91BFDB", "#A6D96A", "#FFD580"])  # Dark Pink, Light Sky Blue, Soft Green, Light Orange

# Plot left hemisphere (lateral view)
plotting.plot_surf_roi(
    surf_mesh=fsaverage["pial_left"],
    roi_map=texture_left,
    hemi="left",
    view="lateral",
    cmap=cmap,
    title="Left Hemisphere (Lateral)"
)

# Plot left hemisphere (medial view)
plotting.plot_surf_roi(
    surf_mesh=fsaverage["pial_left"],
    roi_map=texture_left,
    hemi="left",
    view="medial",
    cmap=cmap,
    title="Left Hemisphere (Medial)"
)

# Plot right hemisphere (lateral view)
plotting.plot_surf_roi(
    surf_mesh=fsaverage["pial_right"],
    roi_map=texture_right,
    hemi="right",
    view="lateral",
    cmap=cmap,
    title="Right Hemisphere (Lateral)"
)

# Plot right hemisphere (medial view)
plotting.plot_surf_roi(
    surf_mesh=fsaverage["pial_right"],
    roi_map=texture_right,
    hemi="right",
    view="medial",
    cmap=cmap,
    title="Right Hemisphere (Medial)"
)

plotting.show()
