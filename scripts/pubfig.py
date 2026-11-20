"""pubfig - tiny helper for journal-ready matplotlib figures.

Copy this file next to your figure script (keeps the script reproducible for
supplementary material), then:

    import pubfig
    spec = pubfig.setup("nature")                      # rcParams for the journal
    fig, axs = plt.subplots(1, 2, figsize=pubfig.figsize("nature", "double", 0.4))
    ...
    pubfig.label_panels(axs, spec)
    pubfig.audit(fig, "nature", "double")              # prints problems, returns list
    pubfig.save(fig, "out/fig1")                       # pdf + svg + 600-dpi tiff + png

CLI:
    python pubfig.py check FILE [--journal nature --width double]   # audit a saved PDF/TIFF/PNG
    python pubfig.py demo OUTDIR                                     # self-test

Only needs matplotlib + numpy (Pillow ships with matplotlib).
Journal numbers are from author guidelines; re-check the target journal before submission.
"""
import os
import re
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

MM = 1 / 25.4  # mm -> inch

# widths in mm; max_h = max figure height mm; fs = body font pt; min_fs = smallest allowed pt
# label = panel label style ('a' lowercase bold, 'A' uppercase bold)
JOURNALS = {
    "nature":   dict(single=89,  onehalf=120, double=183, max_h=170, fs=7, min_fs=5, label="a"),   # Nature, Nat Comms, Sci Rep
    "science":  dict(single=57,  onehalf=121, double=184, max_h=230, fs=7, min_fs=6, label="A"),
    "cell":     dict(single=85,  onehalf=114, double=174, max_h=225, fs=7, min_fs=6, label="A"),   # Cell Press
    "elsevier": dict(single=90,  onehalf=140, double=190, max_h=240, fs=7, min_fs=6, label="a"),
    "acs":      dict(single=82.5, onehalf=120, double=178, max_h=240, fs=7, min_fs=4.5, label="A"),  # JMC, JCIM, ACS Omega...
    "rsc":      dict(single=83,  onehalf=120, double=171, max_h=233, fs=7, min_fs=6, label="a"),
    "wiley":    dict(single=80,  onehalf=120, double=180, max_h=230, fs=7, min_fs=6, label="a"),
    "springer": dict(single=84,  onehalf=129, double=174, max_h=234, fs=8, min_fs=6, label="a"),
    "plos":     dict(single=132, onehalf=150, double=190, max_h=222, fs=8, min_fs=6, label="A"),   # PLOS: 2.63-7.5 in wide
    "frontiers": dict(single=85, onehalf=130, double=180, max_h=230, fs=8, min_fs=6, label="A"),
    "mdpi":     dict(single=85,  onehalf=130, double=170, max_h=230, fs=8, min_fs=6, label="a"),
    "ieee":     dict(single=88.9, onehalf=130, double=181, max_h=230, fs=8, min_fs=6, label="a"),
    "generic":  dict(single=85,  onehalf=130, double=175, max_h=230, fs=7, min_fs=5, label="a"),
}

# Okabe-Ito: distinguishable under common colour-vision deficiencies.
OKABE_ITO = ["#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#56B4E9", "#F0E442", "#000000"]
GREY = "#7F7F7F"                  # controls / vehicle / non-significant
UP, DOWN = "#D55E00", "#0072B2"   # up/down-regulated convention (vermillion / blue)

BAD_CMAPS = {"jet", "rainbow", "hsv", "gist_rainbow", "nipy_spectral", "gist_ncar", "Spectral", "spectral"}


def spec_for(journal):
    key = (journal or "generic").lower()
    if key not in JOURNALS:
        print(f"[pubfig] unknown journal '{journal}', using 'generic'. Known: {', '.join(JOURNALS)}")
        key = "generic"
    return dict(JOURNALS[key], name=key)


def setup(journal="generic", fontsize=None, font="Arial"):
    """Apply journal rcParams. Returns the spec dict (fs, min_fs, widths, label style)."""
    spec = spec_for(journal)
    fs = fontsize or spec["fs"]
    spec["fs"] = fs
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": [font, "Helvetica", "Liberation Sans", "DejaVu Sans"],
        "mathtext.default": "regular",          # math text in the same sans font
        "font.size": fs,
        "axes.titlesize": fs,
        "axes.labelsize": fs,
        "xtick.labelsize": fs - 1,
        "ytick.labelsize": fs - 1,
        "legend.fontsize": fs - 1,
        "legend.title_fontsize": fs - 1,
        "legend.frameon": False,
        "legend.handlelength": 1.5,
        "axes.linewidth": 0.6,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.prop_cycle": mpl.cycler(color=OKABE_ITO),
        "axes.labelpad": 2.5,
        "lines.linewidth": 1.0,
        "lines.markersize": 3.5,
        "patch.linewidth": 0.6,
        "errorbar.capsize": 2,
        "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "xtick.minor.width": 0.4, "ytick.minor.width": 0.4,
        "xtick.major.size": 3, "ytick.major.size": 3,
        "xtick.minor.size": 1.5, "ytick.minor.size": 1.5,
        "xtick.direction": "out", "ytick.direction": "out",
        "xtick.major.pad": 2, "ytick.major.pad": 2,
        "image.cmap": "viridis",
        "pdf.fonttype": 42,                     # TrueType: editable text, no Type 3 (journals reject Type 3)
        "ps.fonttype": 42,
        "svg.fonttype": "none",                 # keep text as text in SVG
        "savefig.dpi": 600,
        "figure.dpi": 150,
        "figure.constrained_layout.use": True,
    })
    import logging  # Arial missing (Linux)? Liberation/DejaVu are used silently
    logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
    return spec


def figsize(journal="generic", width="single", aspect=0.75, height_mm=None):
    """(w, h) in inches. width: 'single'|'onehalf'|'double' or a number in mm.
    aspect = height/width unless height_mm is given. Height is capped at the journal max."""
    spec = spec_for(journal)
    w = float(width) if isinstance(width, (int, float)) else spec[width]
    h = height_mm if height_mm is not None else w * aspect
    return (w * MM, min(h, spec["max_h"]) * MM)


def label_panels(axes, spec=None, labels=None, offset=(-14, 6), size=None, **kw):
    """Bold panel letters at the top-left of each axes. offset is in points;
    make it more negative in x if a letter collides with tick labels."""
    spec = spec or spec_for("generic")
    axes = list(np.ravel(axes))
    if labels is None:
        base = "abcdefghijklmnopqrstuvwxyz" if spec.get("label", "a") == "a" else "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        labels = base[: len(axes)]
    size = size or spec["fs"] + 1
    return [ax.annotate(lab, xy=(0, 1), xycoords="axes fraction", xytext=offset,
                        textcoords="offset points", fontsize=size, fontweight="bold",
                        ha="left", va="bottom", **kw)
            for ax, lab in zip(axes, labels)]


def p_to_stars(p, thresholds=((1e-4, "****"), (1e-3, "***"), (1e-2, "**"), (0.05, "*"))):
    for t, s in thresholds:
        if p < t:
            return s
    return "ns"


def format_p(p):
    """Exact p for text: 'P = 0.012', 'P < 0.0001'. Many Q1 journals prefer exact values."""
    return "P < 0.0001" if p < 1e-4 else (f"P = {p:.3f}" if p >= 1e-3 else f"P = {p:.1e}")


def sig_bracket(ax, x1, x2, y, p=None, text=None, h=None, lw=0.6, fs=None, exact=False):
    """Significance bracket x1..x2 at height y (data coords). Returns top y (for stacking)."""
    if h is None:
        lo, hi = ax.get_ylim()
        h = (hi - lo) * 0.02
    label = text if text is not None else (format_p(p) if exact else p_to_stars(p))
    ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=lw, c="k", clip_on=False)
    ax.text((x1 + x2) / 2, y + h, label, ha="center", va="bottom",
            fontsize=fs or mpl.rcParams["font.size"] - 1)
    return y + h


def _flatten_rgb(path, dpi):
    """Re-save a raster as RGB without alpha (journals reject alpha channels)."""
    from PIL import Image
    im = Image.open(path)
    if im.mode != "RGB":
        rgba = im.convert("RGBA")
        bg = Image.new("RGB", im.size, (255, 255, 255))
        bg.paste(rgba, mask=rgba.split()[-1])
        im = bg
    kw = {"dpi": (dpi, dpi)}
    if path.lower().endswith((".tif", ".tiff")):
        kw["compression"] = "tiff_lzw"
    im.save(path, **kw)


def save(fig, basepath, formats=("pdf", "svg", "tiff", "png"), dpi=600):
    """Save at the exact figure size (no bbox_inches='tight', so the width stays journal-exact).
    dpi: 300 photo/halftone, 600 combination (default), 1000-1200 pure line art."""
    d = os.path.dirname(basepath)
    if d:
        os.makedirs(d, exist_ok=True)
    paths = []
    for ext in formats:
        p = f"{basepath}.{ext}"
        raster = ext in ("tiff", "tif", "png", "jpg")
        fig.savefig(p, dpi=dpi if raster else None, facecolor="white")
        if ext in ("tiff", "tif", "png"):
            _flatten_rgb(p, dpi)
        paths.append(p)
    print("[pubfig] saved:", ", ".join(paths))
    return paths


def audit(fig, journal="generic", width=None, verbose=True):
    """Check a live figure against the journal spec. Returns a list of problem strings."""
    spec = spec_for(journal)
    probs = []
    w_mm, h_mm = (v * 25.4 for v in fig.get_size_inches())
    if width is not None:
        target = float(width) if isinstance(width, (int, float)) else spec[width]
        if abs(w_mm - target) > 1.0:
            probs.append(f"width {w_mm:.1f} mm != {spec['name']} {width} ({target} mm)")
    elif not any(abs(w_mm - spec[k]) <= 1.0 for k in ("single", "onehalf", "double")):
        probs.append(f"width {w_mm:.1f} mm matches no {spec['name']} column width "
                     f"({spec['single']}/{spec['onehalf']}/{spec['double']} mm)")
    if h_mm > spec["max_h"] + 0.5:
        probs.append(f"height {h_mm:.1f} mm exceeds max {spec['max_h']} mm")

    fig.canvas.draw()
    small = {f"'{t.get_text()[:20]}' {t.get_fontsize():.1f}pt"
             for t in fig.findobj(mpl.text.Text)
             if t.get_visible() and t.get_text().strip() and t.get_fontsize() < spec["min_fs"] - 1e-6}
    if small:
        probs.append(f"text below {spec['min_fs']} pt: " + "; ".join(sorted(small)[:6]))
    # mathtext sub/superscripts ($_{50}$, $^{-1}$) render at 70 % of the label size
    scripts = {f"'{t.get_text()[:20]}' -> {0.7 * t.get_fontsize():.1f}pt"
               for t in fig.findobj(mpl.text.Text)
               if t.get_visible() and re.search(r"\$[^$]*[_^]", t.get_text())
               and 0.7 * t.get_fontsize() < spec["min_fs"] - 1e-6}
    notes = []
    if scripts:
        notes.append(f"sub/superscripts below {spec['min_fs']} pt (usually tolerated; for strict checks raise "
                     f"those labels to {spec['min_fs'] / 0.7:.1f} pt): " + "; ".join(sorted(scripts)[:4]))
    thin = [l for l in fig.findobj(mpl.lines.Line2D)
            if l.get_visible() and 0 < l.get_linewidth() < 0.25 and l.get_linestyle() not in ("None", "")]
    if thin:
        probs.append(f"{len(thin)} line(s) thinner than 0.25 pt (vanish in print)")
    for sm in fig.findobj(lambda a: isinstance(a, mpl.cm.ScalarMappable)):
        cm = getattr(sm, "cmap", None)
        if cm is not None and cm.name.removesuffix("_r") in BAD_CMAPS:
            probs.append(f"rainbow colormap '{cm.name}' (not perceptually uniform / CVD-safe); use viridis, cividis or RdBu_r")
    unlabeled = [ax for ax in fig.axes if ax.axison and ax.name == "rectilinear"
                 and ax.get_label() != "<colorbar>" and not (ax.get_xlabel() or ax.get_ylabel())]
    if unlabeled:
        probs.append(f"{len(unlabeled)} axes without any axis label")
    if verbose:
        if probs:
            print("[pubfig] AUDIT ISSUES:")
            for p in probs:
                print("   -", p)
        else:
            print(f"[pubfig] audit OK ({w_mm:.1f} x {h_mm:.1f} mm, {spec['name']})")
        for n in notes:
            print("[pubfig] note:", n)
    return probs


def check_file(path, journal="generic", width=None):
    """Audit a saved file: physical size, dpi + colour mode (raster), Type 3 fonts (PDF)."""
    spec = spec_for(journal)
    probs, info = [], {}
    if path.lower().endswith(".pdf"):
        raw = open(path, "rb").read()
        m = re.search(rb"/MediaBox\s*\[\s*([\d.\-]+)\s+([\d.\-]+)\s+([\d.\-]+)\s+([\d.\-]+)", raw)
        if m:
            x0, y0, x1, y1 = map(float, m.groups())
            info["w_mm"], info["h_mm"] = (x1 - x0) / 72 * 25.4, (y1 - y0) / 72 * 25.4
        info["type3_fonts"] = b"/Type3" in raw
        if info["type3_fonts"]:
            probs.append("PDF contains Type 3 fonts (set pdf.fonttype=42)")
    else:
        from PIL import Image
        im = Image.open(path)
        dpi = float(im.info.get("dpi", (72, 72))[0])
        info.update(px=im.size, dpi=round(dpi), mode=im.mode,
                    w_mm=im.size[0] / dpi * 25.4, h_mm=im.size[1] / dpi * 25.4)
        if dpi < 300:
            probs.append(f"raster dpi {dpi:.0f} < 300")
        if im.mode not in ("RGB", "CMYK", "L"):
            probs.append(f"raster mode {im.mode} (flatten to RGB)")
    if "w_mm" in info:
        if width is not None:
            target = spec[width] if width in spec else float(width)
            if abs(info["w_mm"] - target) > 1.0:
                probs.append(f"width {info['w_mm']:.1f} mm != {target} mm")
        if info["h_mm"] > spec["max_h"] + 0.5:
            probs.append(f"height {info['h_mm']:.1f} mm > max {spec['max_h']} mm")
    return info, probs


def _demo(outdir):
    spec = setup("nature")
    rng = np.random.default_rng(0)
    fig, axs = plt.subplots(1, 2, figsize=figsize("nature", "double", 0.38))
    groups = ["Vehicle", "Drug A", "Drug B"]
    data = [rng.normal(m, 8, 6) for m in (100, 70, 45)]
    axs[0].boxplot(data, widths=0.5, showfliers=False, medianprops=dict(color="k"))
    for i, d in enumerate(data, 1):
        axs[0].scatter(np.full(d.size, i) + rng.uniform(-.12, .12, d.size), d, s=8,
                       color=[GREY, OKABE_ITO[0], OKABE_ITO[3]][i - 1], zorder=3, lw=0)
    axs[0].set_xticks([1, 2, 3], groups)
    axs[0].set_ylabel("Cell viability (%)")
    sig_bracket(axs[0], 1, 3, 118, p=1e-5)
    x = np.linspace(0, 100, 500)
    for i in range(3):
        axs[1].plot(x, 0.15 + 0.1 * i + 0.02 * np.sin(x / 5 + i) + 0.01 * rng.standard_normal(500),
                    lw=0.8, label=f"Complex {i + 1}")
    axs[1].set_xlabel("Time (ns)")
    axs[1].set_ylabel("RMSD (nm)")
    axs[1].legend(loc="upper left", ncol=3)
    label_panels(axs, spec)
    probs = audit(fig, "nature", "double")
    paths = save(fig, os.path.join(outdir, "demo"))
    for p in paths:
        if p.endswith((".pdf", ".tiff")):
            info, fp = check_file(p, "nature", "double")
            print("  ", os.path.basename(p), info, fp or "OK")
            assert not fp, fp
    assert not probs, probs
    assert p_to_stars(0.03) == "*" and p_to_stars(0.2) == "ns" and p_to_stars(5e-5) == "****"
    print("[pubfig] demo self-check passed")


if __name__ == "__main__":
    a = sys.argv[1:]
    opt = lambda k: a[a.index(k) + 1] if k in a else None  # noqa: E731
    if a[:1] == ["demo"]:
        _demo(a[1] if len(a) > 1 else ".")
    elif a[:1] == ["check"] and len(a) > 1:
        j, w = opt("--journal") or "generic", opt("--width")
        files = [x for i, x in enumerate(a[1:], 1) if not x.startswith("--") and a[i - 1] not in ("--journal", "--width")]
        bad = 0
        for f in files:
            info, probs = check_file(f, j, w)
            print(f, {k: (round(v, 1) if isinstance(v, float) else v) for k, v in info.items()})
            print("   ", "OK" if not probs else "; ".join(probs))
            bad += bool(probs)
        sys.exit(1 if bad else 0)
    else:
        print(__doc__)
