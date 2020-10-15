import cv2
import plotly.express as px
import os
from tqdm import tqdm
import tensorflow_hub as hub
import tensorflow as tf
from PIL import Image

# For drawing onto the image.
import numpy as np
from PIL import Image
from PIL import ImageColor
from PIL import ImageDraw
from PIL import ImageFont
from PIL import ImageOps


def draw_bounding_box_on_image(
    image, ymin, xmin, ymax, xmax, color, font, thickness=4, display_str_list=()
):
    """Adds a bounding box to an image."""
    draw = ImageDraw.Draw(image)
    im_width, im_height = image.size
    (left, right, top, bottom) = (
        xmin * im_width,
        xmax * im_width,
        ymin * im_height,
        ymax * im_height,
    )
    draw.line(
        [(left, top), (left, bottom), (right, bottom), (right, top), (left, top)],
        width=thickness,
        fill=color,
    )

    # If the total height of the display strings added to the top of the bounding
    # box exceeds the top of the image, stack the strings below the bounding box
    # instead of above.
    display_str_heights = [font.getsize(ds)[1] for ds in display_str_list]
    # Each display_str has a top and bottom margin of 0.05x.
    total_display_str_height = (1 + 2 * 0.05) * sum(display_str_heights)

    if top > total_display_str_height:
        text_bottom = top
    else:
        text_bottom = top + total_display_str_height
    # Reverse list and print from bottom to top.
    for display_str in display_str_list[::-1]:
        text_width, text_height = font.getsize(display_str)
        margin = np.ceil(0.05 * text_height)
        draw.rectangle(
            [
                (left, text_bottom - text_height - 2 * margin),
                (left + text_width, text_bottom),
            ],
            fill=color,
        )
        draw.text(
            (left + margin, text_bottom - text_height - margin),
            display_str,
            fill="black",
            font=font,
        )
        text_bottom -= text_height - 2 * margin


def draw_boxes(
    image, boxes, class_ids, class_names, scores, font, max_boxes=10, min_score=0.1
):
    """Overlay labeled boxes on an image with formatted scores and label names."""
    colors = list(ImageColor.colormap.values())

    for i in range(min(boxes.shape[0], max_boxes)):
        if scores[i] >= min_score:
            ymin, xmin, ymax, xmax = tuple(boxes[i])
            display_str = "{}: {}%".format(
                class_names[i].decode("ascii"), int(100 * scores[i])
            )
            color = colors[class_ids[i] % len(colors)]
            image_pil = Image.fromarray(np.uint8(image)).convert("RGB")
            draw_bounding_box_on_image(
                image_pil,
                ymin,
                xmin,
                ymax,
                xmax,
                color,
                font,
                display_str_list=[display_str],
            )
            np.copyto(image, np.array(image_pil))
    return np.array(image_pil)


def fast_draw_boxes(
    image, boxes, class_ids, class_names, scores, font, max_boxes=10, min_score=0.1
):
    """Overlay labeled boxes on an image with formatted scores and label names."""
    colors = list(ImageColor.colormap.values())
    image_pil = Image.fromarray(np.uint8(image)).convert("RGB")

    for i in range(min(boxes.shape[0], max_boxes)):
        if scores[i] >= min_score:
            ymin, xmin, ymax, xmax = tuple(boxes[i])
            display_str = "{}: {}%".format(
                class_names[i].decode("ascii"), int(100 * scores[i])
            )
            color = colors[class_ids[i] % len(colors)]
            draw_bounding_box_on_image(
                image_pil,
                ymin,
                xmin,
                ymax,
                xmax,
                color,
                font,
                display_str_list=[display_str],
            )
    #       np.copyto(image, np.array(image_pil))
    return np.array(image_pil)


module_handle = "https://tfhub.dev/google/openimages_v4/ssd/mobilenet_v2/1"
detector = hub.load(module_handle).signatures["default"]
codec = "XVID"

try:
    font = ImageFont.truetype(
        "/usr/share/fonts/truetype/liberation/LiberationSansNarrow-Regular.ttf", 25
    )
except IOError:
    print("Font not found, using default font.")
    font = ImageFont.load_default()

i = 2
# Define the codec and create VideoWriter object
VIDEO_PATH = f"./data/scene-{i}.mov"
VIDEO_OUT = f"./data/processed/scene_{i}.mov"

frames = []

cap = cv2.VideoCapture(VIDEO_PATH)
ret = True

while ret:
    ret, frame = cap.read()
    if ret:
        frame = cv2.rotate(frame, cv2.ROTATE_90_COUNTERCLOCKWISE)
        frame = cv2.cvtColor(frame, cv2.COLOR_BGR2RGB)
        frames.append(frame)

processed_frames = []
sample_rate = 3

for i, img in enumerate(tqdm(frames)):
    if i % sample_rate == 0:
        resized = cv2.resize(img, (512, 512))
        image_tensor = tf.image.convert_image_dtype(resized, tf.float32)[
            tf.newaxis, ...
        ]

        result = detector(image_tensor)

    image_with_boxes = fast_draw_boxes(
        img.copy(),
        result["detection_boxes"].numpy(),
        result["detection_class_labels"].numpy(),
        result["detection_class_entities"].numpy(),
        result["detection_scores"].numpy(),
        font=font,
    )

    processed_frames.append(image_with_boxes)

fourcc = cv2.VideoWriter_fourcc(*codec)
out = cv2.VideoWriter(VIDEO_OUT, fourcc, 30, (1280, 720))

for frame in tqdm(processed_frames):
    frame = cv2.cvtColor(frame, cv2.COLOR_RGB2BGR)
    out.write(frame)

out.release()
