from joblib import Memory
from utils.trainable_segmentation import multiscale_basic_features
from skimage import io as skio


memory = Memory("./joblib_cache", bytes_limit=3000000000, verbose=3)

compute_features = memory.cache(multiscale_basic_features)

DEFAULT_STROKE_WIDTH = 3  # gives line width of 2^3 = 8

DEFAULT_IMAGE_PATH = "assets/images/segmentation_img.jpg"

SEG_FEATURE_TYPES = ["intensity", "edges", "texture"]

# the number of different classes for labels
NUM_LABEL_CLASSES = 5
DEFAULT_LABEL_CLASS = 0
class_label_colormap = ["#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2"]
class_labels = list(range(NUM_LABEL_CLASSES))
# we can't have less colors than classes
assert NUM_LABEL_CLASSES <= len(class_label_colormap)

# Font and background colors associated with each theme
text_color = {"dark": "#95969A", "light": "#595959"}
card_color = {"dark": "#2D3038", "light": "#FFFFFF"}

img = skio.imread(DEFAULT_IMAGE_PATH)
