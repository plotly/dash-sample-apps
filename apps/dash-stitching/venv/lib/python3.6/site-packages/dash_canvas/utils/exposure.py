import numpy as np

def contrast_adjust(img, intensity):
    print("contrast adjust", intensity)
    contrast_range = 255 ** (2 - 2*intensity)
    img_f = img.astype(np.float)
    img_f = np.clip((img_f - 127) * 255. / contrast_range + 127, 0, 255)
    return img_f.astype(np.uint8)


def brightness_adjust(img, intensity):
    offset = int((intensity - 0.5) * 255)
    if intensity < 0.5:
        img = np.clip(img, np.abs(offset), 255)
        return (img + offset).astype(np.uint8)
    elif intensity > 0.5:
        img = np.clip(img, None, 255 - offset)
        return (img + offset).astype(np.uint8)
    else:
        return img

