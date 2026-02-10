"""
Cocofluo pipeline configuration.
"""

N_CHANNELS = 13  # Number of channels in the image
DAPI_CH = 5  # Structural channel index (1-based)

CHANNEL_NAMES = {
    1: "adcyap1",
    2: "bdhf",
    3: "casr",
    4: "cck",
    5: "dapi",
    6: "dbh",
    7: "fos",
    8: "gaba",
    9: "glp1r",
    10: "htr2c",
    11: "htr3a",
    12: "tac1",
    13: "vglut2",
}

DAPI_ENHANCEMENT = True  # Apply contrast stretch (1–99 percentile normalization)

USE_GPU = True
FLOW_THRESHOLD = 0.4
CELLPROB_THRESHOLD = 0.0

MARKER_EXTRACT_PERCENTILE = 99.8  # Top (100-P)% pixels = positive;

MARKER_EXPRESS_THRESHOLD = 0.08  # ratio of marker-positive area / cell area

COLOCALIZATION_GROUPS = [
    (3, 7),
    (7, 9),
    (3, 7, 9),
]

CELL_SEG_SAMPLING = True  # Save cell segmentation QC images
MARKER_EXTRACT_SAMPLING = True  # Save marker extraction QC images
MARKER_EXPRESS_SAMPLING = True  # Save express-cell highlight images

INPUT_DIR = "source"
SAMPLE_DIR = "sampling"
OUTPUT_DIR = "results"
