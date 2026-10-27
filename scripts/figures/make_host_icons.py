import json, os
from PIL import Image
import numpy as np

CROPS = r"D:\Claude_workzone\test\sil_crops"
OUT = r"D:\Claude_workzone\icons"
os.makedirs(OUT, exist_ok=True)

# contact-sheet index -> species (visually checked)
MAP = {
    "Ovis_aries": 8, "Pantholops_hodgsonii": 13, "Bos_taurus": 0, "Bos_grunniens": 5,
    "Capra_hircus": 7, "Bubalus_bubalis": 3, "Bos_indicus": 2, "Procapra_picticaudata": 16,
    "Bos_frontalis": 4,
    "Hydropotes_inermis": 9, "Alces_alces": 14, "Cervus_elaphus": 12,
    "Rangifer_tarandus": 11, "Alces_americanus": 15,
    "Giraffa_camelopardalis": 19,
    "Camelus_bactrianus": 18, "Camelus_dromedarius": 17,
}

for name, idx in MAP.items():
    rgb = np.asarray(Image.open(os.path.join(CROPS, f"{idx:02d}.png")).convert("L"), dtype=float)
    # white/off-white background -> transparent; black silhouette -> opaque
    alpha = np.clip((235.0 - rgb) / (235.0 - 40.0), 0, 1) * 255
    rgba = np.zeros((*rgb.shape, 4), dtype=np.uint8)
    rgba[..., :3] = 0x3B, 0x42, 0x52  # infographic dark slate
    rgba[..., 3] = alpha.astype(np.uint8)
    img = Image.fromarray(rgba, "RGBA")
    bbox = img.getbbox()
    img.crop(bbox).save(os.path.join(OUT, f"{name}.png"))
    print(name, img.crop(bbox).size)
