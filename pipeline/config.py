"""
Cocofluo pipeline configuration.
"""

N_CHANNELS = 4  # Number of channels in the image
DAPI_CH = 1  # Structural channel index (1-based)

CHANNEL_NAMES = {
    1: "dapi",
    2: "cas",
    3: "fos",
    4: "glp1r",
    5: "",
    6: "",
    7: "",
    8: "",
    9: "",
    10: "",
    11: "",
    12: "",
}

DAPI_ENHANCEMENT = True  # Apply contrast stretch (1–99 percentile normalization)

USE_GPU = True
FLOW_THRESHOLD = 0.4
CELLPROB_THRESHOLD = 0.0

MARKER_EXTRACT_PERCENTILE = 99.8  # Top (100-P)% pixels = positive;

MARKER_EXPRESS_THRESHOLD = 0.07  # ratio of marker-positive area / cell area

COLOCALIZATION_GROUPS = [
    (2, 4),
    (3, 4),
    (2, 3, 4),
]

CELL_SEG_SAMPLING = True  # Save cell segmentation QC images
MARKER_EXTRACT_SAMPLING = True  # Save marker extraction QC images
MARKER_EXPRESS_SAMPLING = True  # Save express-cell highlight images

INPUT_DIR = "source"
SAMPLE_DIR = "sampling"
OUTPUT_DIR = "results"
