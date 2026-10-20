#!/usr/bin/env python3
"""Cut the README figures out of the unified-annotation slide deck.

The deck (figs/unified_annotation_scheme.odp, first built by
scripts/make_unified_slides.js and then edited by hand) has 5 slides:
slide 1 shows the tiers, slides 2-5 the four special cases. For the README
this writes two figures:

  figs/unified_annotation_tiers.png  slide 1
  figs/unified_annotation_cases.png  slides 2-5 as a 2 x 2 grid

Export the deck to PDF first (PowerPoint / Impress "Export as PDF", or
`soffice --headless --convert-to pdf figs/unified_annotation_scheme.odp`),
then:

  scripts/make_unified_readme_figures.py unified_annotation_scheme.pdf [-o figs]

Needs `pdftoppm` (poppler) and Pillow. Dev-only; not used by the pipeline.
"""
import argparse
import os
import subprocess
import tempfile

from PIL import Image, ImageChops, ImageDraw

DPI = 150
SLIDE_NUMBER_BOX = (0.93, 0.92)  # bottom-right corner (fraction of w, h) holding the slide number
PAD = 30  # px of white kept around the cropped content


def render(pdf, dpi, outdir):
    subprocess.run(["pdftoppm", "-r", str(dpi), "-png", pdf, os.path.join(outdir, "p")], check=True)
    pages = sorted(f for f in os.listdir(outdir) if f.startswith("p-"))
    return [Image.open(os.path.join(outdir, f)).convert("RGB") for f in pages]


def blank_slide_number(im):
    w, h = im.size
    ImageDraw.Draw(im).rectangle([int(w * SLIDE_NUMBER_BOX[0]), int(h * SLIDE_NUMBER_BOX[1]), w, h], fill="white")
    return im


def content_box(im):
    bg = Image.new("RGB", im.size, "white")
    return ImageChops.difference(im, bg).convert("L").point(lambda v: 255 if v > 8 else 0).getbbox()


def crop(im, box):
    x0, y0, x1, y1 = box
    w, h = im.size
    return im.crop((max(0, x0 - PAD), max(0, y0 - PAD), min(w, x1 + PAD), min(h, y1 + PAD)))


def union(boxes):
    return (min(b[0] for b in boxes), min(b[1] for b in boxes),
            max(b[2] for b in boxes), max(b[3] for b in boxes))


def save(im, path):
    tmp = path + ".tmp.png"
    im.save(tmp, optimize=True)
    os.replace(tmp, path)
    print(f"wrote {path} ({im.size[0]} x {im.size[1]} px)")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("pdf", help="PDF export of the 5-slide deck")
    ap.add_argument("-o", "--outdir", default="figs")
    ap.add_argument("--dpi", type=int, default=DPI)
    a = ap.parse_args()

    with tempfile.TemporaryDirectory() as td:
        slides = [blank_slide_number(s) for s in render(a.pdf, a.dpi, td)]
    if len(slides) != 5:
        raise SystemExit(f"expected 5 slides, got {len(slides)}")

    save(crop(slides[0], content_box(slides[0])), os.path.join(a.outdir, "unified_annotation_tiers.png"))

    # One shared crop box for the four case slides keeps the panels the same
    # size and their tracks aligned across the grid.
    cases = slides[1:]
    box = union([content_box(s) for s in cases])
    panels = [crop(s, box) for s in cases]
    pw, ph = panels[0].size
    gap = 4
    grid = Image.new("RGB", (2 * pw + gap, 2 * ph + gap), (200, 200, 200))
    for i, p in enumerate(panels):
        grid.paste(p, ((i % 2) * (pw + gap), (i // 2) * (ph + gap)))
    save(grid, os.path.join(a.outdir, "unified_annotation_cases.png"))


if __name__ == "__main__":
    main()
