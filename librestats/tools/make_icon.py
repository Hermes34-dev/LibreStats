"""Make the Extension Manager icon (src/icons/librestats.png) from the logo.

Development tool only; needs Pillow (not shipped):
    python tools/make_icon.py ../LibreStats_LOGO.png src/icons/librestats.png

Takes the blue document shape (the logo without its text), makes the white
background outside it transparent, centres it on a square and scales it down.
"""

import sys

from PIL import Image, ImageDraw

SIZE = 42  # size the Extension Manager shows


def main(src, dst):
    logo = Image.open(src).convert("RGBA")
    w, h = logo.size
    # The document shape: everything above the "LibreStats" text (top ~79 %).
    doc = logo.crop((0, 0, w, int(h * 0.79)))
    bbox = doc.convert("L").point(lambda v: 255 if v < 200 else 0).getbbox()
    doc = doc.crop(bbox)
    # Clear the white surround: flood-fill from each corner. The white chart inside the
    # shape is enclosed by blue, so it stays.
    for corner in ((0, 0), (doc.width - 1, 0), (0, doc.height - 1),
                   (doc.width - 1, doc.height - 1)):
        if doc.getpixel(corner)[3] and min(doc.getpixel(corner)[:3]) > 200:
            ImageDraw.floodfill(doc, corner, (255, 255, 255, 0), thresh=60)
    side = max(doc.size)
    square = Image.new("RGBA", (side, side), (255, 255, 255, 0))
    square.paste(doc, ((side - doc.width) // 2, (side - doc.height) // 2))
    square.resize((SIZE, SIZE), Image.LANCZOS).save(dst, optimize=True)


if __name__ == "__main__":
    main(*sys.argv[1:3])
