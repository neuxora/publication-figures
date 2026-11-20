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

Chemical structures (needs rdkit): pubfig.draw_structures(fig, smiles, labels, core=scaffold)

Only needs matplotlib + numpy (Pillow ships with matplotlib); rdkit only for structures.
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


# ---- chemical structures (needs rdkit: pip install rdkit) --------------------------
ACS1996_PT = 10  # ACS Document 1996: bond 14.4 pt (0.508 cm), line 0.6 pt, Arial 10 pt labels


def _highlights(mol, highlight):
    """{SMARTS or atom-index tuple: colour} (or a list -> Okabe-Ito colours) -> RDKit highlight args."""
    from rdkit import Chem
    if not highlight:
        return {}
    if not isinstance(highlight, dict):
        highlight = dict(zip(highlight, [OKABE_ITO[i] for i in (1, 2, 4, 3, 5, 0)]))
    atoms, bonds = {}, {}
    for key, colour in highlight.items():
        rgb = tuple(mpl.colors.to_rgb(colour))
        if isinstance(key, str):
            patt = Chem.MolFromSmarts(key)
            if patt is None:
                raise ValueError(f"invalid SMARTS: {key}")
            matches = mol.GetSubstructMatches(patt)
        else:
            matches = [tuple(key)]
        for match in matches:
            atoms.update({a: rgb for a in match})
            bonds.update({b.GetIdx(): rgb for b in mol.GetBonds()
                          if b.GetBeginAtomIdx() in match and b.GetEndAtomIdx() in match})
    return dict(highlightAtoms=list(atoms), highlightBonds=list(bonds),
                highlightAtomColors=atoms, highlightBondColors=bonds)


def mol_svg(mol, highlight=None, highlight_style="halo"):
    """SVG of a molecule (RDKit Mol or SMILES) in ACS Document 1996 style; 1 SVG unit = 1 pt.

    highlight: {SMARTS or atom indices: colour} or a list of SMARTS (Okabe-Ito colours), e.g.
    {"C(F)(F)F": "#E69F00"}. highlight_style "halo" = light tint behind the group (structure stays
    black, survives greyscale print); "color" = the group's bonds and labels drawn in colour.
    """
    from rdkit import Chem
    from rdkit.Chem import Draw, rdDepictor
    from rdkit.Chem.Draw import rdMolDraw2D
    if isinstance(mol, str):
        smi, mol = mol, Chem.MolFromSmiles(mol)
        if mol is None:
            raise ValueError(f"invalid SMILES: {smi}")
    if not mol.GetNumConformers():
        rdDepictor.SetPreferCoordGen(False)  # CoordGen + ACS mode draws bonds through charged labels
        rdDepictor.Compute2DCoords(mol)
    mol = Chem.Mol(mol)
    for a in mol.GetAtoms():  # [*:1] in a scaffold SMILES -> R1 with subscript
        if a.GetAtomicNum() == 0 and a.GetAtomMapNum():
            a.SetProp("atomLabel", f"R<sub>{a.GetAtomMapNum()}</sub>")
            a.SetAtomMapNum(0)
    d = rdMolDraw2D.MolDraw2DSVG(-1, -1)
    o = d.drawOptions()
    o.fontFile = mpl.font_manager.findfont("Arial")  # falls back to a metric-compatible sans
    rdMolDraw2D.SetACS1996Mode(o, Draw.MeanBondLength(mol) or 1.5)
    o.clearBackground = False
    hl = _highlights(mol, highlight)
    if hl and highlight_style not in ("halo", "color"):
        raise ValueError("highlight_style must be 'halo' or 'color'")
    if hl and highlight_style == "color":  # the group's bonds and labels in colour
        o.circleAtoms = False
        o.continuousHighlight = False
        d.DrawMolecule(mol, **hl)
    else:
        d.DrawMolecule(mol)
    d.FinishDrawing()
    svg = d.GetDrawingText()
    if hl and highlight_style == "halo":  # marker-pen band under the black structure
        pos = [d.GetDrawCoords(i) for i in range(mol.GetNumAtoms())]
        lens = [(pos[b.GetBeginAtomIdx()] - pos[b.GetEndAtomIdx()]).Length() for b in mol.GetBonds()]
        w = 0.6 * (float(np.median(lens)) if lens else 14.4)
        hexc = lambda c: mpl.colors.to_hex(tuple(0.45 * v + 0.55 for v in c))  # 55 % tint
        els = []
        for b, c in hl["highlightBondColors"].items():
            p, q = (pos[i] for i in (mol.GetBondWithIdx(b).GetBeginAtomIdx(), mol.GetBondWithIdx(b).GetEndAtomIdx()))
            els.append(f"<path class='halo' d='M {p.x:.2f},{p.y:.2f} L {q.x:.2f},{q.y:.2f}' "
                       f"style='fill:none;stroke:{hexc(c)};stroke-width:{w:.2f}px;stroke-linecap:round'/>")
        for a, c in hl["highlightAtomColors"].items():
            els.append(f"<ellipse class='halo' cx='{pos[a].x:.2f}' cy='{pos[a].y:.2f}' rx='{w / 2:.2f}' "
                       f"ry='{w / 2:.2f}' style='fill:{hexc(c)};stroke:none'/>")
        svg = svg.replace("<!-- END OF HEADER -->", "<!-- END OF HEADER -->\n" + "\n".join(els), 1)
    if _label_clashes(svg):
        import warnings
        warnings.warn(f"[pubfig] bond drawn through an atom label in {Chem.MolToSmiles(mol)}; "
                      "check the render or supply 2D coordinates (e.g. a MOL file from ChemDraw)")
    return svg


def _label_clashes(svg, inset=0.8):
    """Count bond end points inside a label glyph's box (guards against RDKit drawing bugs)."""
    boxes, ends = [], []
    for el in re.findall(r"<path[^>]*>", svg, re.S):
        cls = re.search(r"class='([^']*)'", el)
        cls = cls.group(1) if cls else ""
        nums = [float(v) for v in re.findall(r"-?\d+\.?\d*", re.search(r" d='([^']*)'", el).group(1))]
        pts = list(zip(nums[0::2], nums[1::2]))
        if re.fullmatch(r"atom-\d+", cls):  # one path per glyph
            xs, ys = zip(*pts)
            boxes.append((min(xs), min(ys), max(xs), max(ys)))
        elif cls.startswith("bond-") and "fill:none" in el:  # stroked bond lines, not halos/wedges
            ends += [pts[0], pts[-1]]
    return sum(b[0] + inset < x < b[2] - inset and b[1] + inset < y < b[3] - inset
               for x, y in ends for b in boxes)


def _svg_path(d):
    """RDKit path data (absolute M/L/Q/Z) -> matplotlib Path."""
    from matplotlib.path import Path
    code = {"M": Path.MOVETO, "L": Path.LINETO, "Q": Path.CURVE3}
    verts, codes, nums, cmd = [], [], [], "M"
    for tok in re.findall(r"[MLQZ]|-?[\d.]+(?:e-?\d+)?", d):
        if tok in code:
            cmd = tok
        elif tok == "Z":
            if verts:
                verts.append(verts[-1]); codes.append(Path.CLOSEPOLY)
        else:
            nums.append(float(tok))
            if len(nums) == 2:
                verts.append(tuple(nums)); codes.append(code[cmd]); nums = []
                if cmd == "M":
                    cmd = "L"  # implicit lineto after moveto
    return Path(verts, codes)


def draw_structure(fig, mol, x_mm, y_mm, label=None, scale=1.0, label_size=None,
                   highlight=None, highlight_style="halo"):
    """Draw a vector structure centred at (x_mm, y_mm), measured from the figure's bottom-left.

    scale=1 is ACS 1996 at 100 % (bond 14.4 pt, atom labels 10 pt). Below 0.6 the atom
    labels drop under 6 pt; audit() flags it. `label` (e.g. "4a") is drawn bold beneath.
    highlight / highlight_style: see mol_svg. Returns the axes (printed size: ax.structure_mm).
    """
    from matplotlib.patches import Ellipse, PathPatch
    svg = (mol if isinstance(mol, str) and mol.lstrip().startswith("<?xml")
           else mol_svg(mol, highlight, highlight_style))
    W, H = _svg_size(svg)
    if not any(a.get_subplotspec() for a in fig.axes):
        fig.set_layout_engine("none")  # structures-only figure: nothing for constrained layout to do
    fw, fh = fig.get_size_inches() * 72
    w, h = W * scale, H * scale
    x0, y0 = x_mm / 25.4 * 72 - w / 2, y_mm / 25.4 * 72 - h / 2
    ax = fig.add_axes([x0 / fw, y0 / fh, w / fw, h / fh], label=f"<structure {scale:g}>")
    ax.set_axis_off()
    ax.set_in_layout(False)
    ax.set_xlim(0, W); ax.set_ylim(H, 0)
    ax.highlighted = "<ellipse" in svg or bool(re.search(r"(?:fill|stroke)[:=]'?#(?!000000|FFFFFF)", svg))
    for el in re.findall(r"<(?:path|ellipse)[^>]*>", svg, re.S):
        style = re.search(r"style='([^']*)'", el)
        st = dict(kv.split(":", 1) for kv in (style.group(1) if style else "").split(";") if ":" in kv)
        fill = re.search(r" fill='([^']*)'", el)
        fill = fill.group(1) if fill else st.get("fill", "none")
        if el.startswith("<ellipse"):  # highlight halo on an atom
            g = {k: float(re.search(rf" {k}='([\d.\-]+)'", el).group(1)) for k in ("cx", "cy", "rx", "ry")}
            ax.add_patch(Ellipse((g["cx"], g["cy"]), 2 * g["rx"], 2 * g["ry"], fc=fill, ec="none", lw=0))
            continue
        path = _svg_path(re.search(r" d='([^']*)'", el).group(1))
        if fill != "none":
            ax.add_patch(PathPatch(path, fc=fill, ec="none", lw=0))
        if st.get("stroke", "none") != "none":
            cap = {"round": "round", "square": "projecting"}.get(st.get("stroke-linecap"), "butt")
            ax.add_patch(PathPatch(path, fill=False, ec=st["stroke"], capstyle=cap, joinstyle="miter",
                                   lw=float(st.get("stroke-width", "0.6").rstrip("px")) * scale))
    if label:
        ax.annotate(label, (W / 2, H), xytext=(0, -3), textcoords="offset points", ha="center",
                    va="top", fontweight="bold", annotation_clip=False,
                    size=label_size or max(ACS1996_PT * scale, mpl.rcParams["font.size"]))
    ax.structure_mm = (w / 72 * 25.4, h / 72 * 25.4)
    return ax


def _svg_size(svg):
    return tuple(map(float, re.search(r"viewBox='0 0 ([\d.]+) ([\d.]+)'", svg).groups()))


def draw_structures(fig, mols, labels=None, ncols=None, rect_mm=None, gap_mm=4, max_scale=1.0, core=None,
                    highlight=None, highlight_style="halo"):
    """Grid of structures at ONE common scale (equal bond lengths, as journals expect).

    rect_mm = (x, y, w, h) area in mm from the figure's bottom-left (default: whole figure),
    or an empty panel axes from subplots/subplot_mosaic to fill (multi-panel figures).
    Scale is the largest <= max_scale that fits every cell; labels share a baseline per row.
    core = SMILES/SMARTS of the shared scaffold: every compound containing it is drawn with
    the scaffold in the same orientation (series figures). highlight = {SMARTS: colour} applied
    to every compound that contains the group (see mol_svg). Returns (list_of_axes, scale).
    """
    if core is not None:
        from rdkit import Chem
        from rdkit.Chem import rdDepictor
        rdDepictor.SetPreferCoordGen(False)
        core = Chem.MolFromSmiles(core) or Chem.MolFromSmarts(core) if isinstance(core, str) else core
        rdDepictor.Compute2DCoords(core)
        mols = [Chem.MolFromSmiles(m) if isinstance(m, str) else Chem.Mol(m) for m in mols]
        for m in mols:
            if m.HasSubstructMatch(core):
                rdDepictor.GenerateDepictionMatching2DStructure(m, core)
    svgs = [mol_svg(m, highlight, highlight_style) for m in mols]
    fw, fh = (v * 25.4 for v in fig.get_size_inches())
    if hasattr(rect_mm, "get_position"):  # an empty panel axes: fill it (after layout settles)
        fig.canvas.draw()
        p = rect_mm.get_position()
        rect_mm.set_axis_off()
        rect_mm = (p.x0 * fw, p.y0 * fh, p.width * fw, p.height * fh)
    x, y, w, h = rect_mm or (0, 0, fw, fh)
    n = len(svgs)
    ncols = ncols or n
    nrows = -(-n // ncols)
    lab_mm = 6 if labels else 0
    cw = (w - gap_mm * (ncols - 1)) / ncols
    ch = (h - gap_mm * (nrows - 1)) / nrows - lab_mm
    pt_mm = 25.4 / 72
    scale = min([max_scale] + [cw / (W * pt_mm) for W, _ in map(_svg_size, svgs)]
                + [ch / (H * pt_mm) for _, H in map(_svg_size, svgs)])
    lsize = max(ACS1996_PT * scale, mpl.rcParams["font.size"])
    axes = []
    for i, svg in enumerate(svgs):
        r, c = divmod(i, ncols)
        cx = x + c * (cw + gap_mm) + cw / 2
        top = y + h - r * (ch + lab_mm + gap_mm)
        axes.append(draw_structure(fig, svg, cx, top - ch / 2, scale=scale))
        if labels:
            fig.text(cx / fw, (top - ch - 1) / fh, labels[i], ha="center", va="top",
                     fontweight="bold", size=lsize)
    return axes, scale


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
    for ax in fig.axes:
        m = re.fullmatch(r"<structure ([\d.]+)>", ax.get_label())
        if not m:
            continue
        if ACS1996_PT * float(m.group(1)) < spec["min_fs"] - 1e-6:
            probs.append(f"structures at scale {float(m.group(1)):.2f}: atom labels {ACS1996_PT * float(m.group(1)):.1f} pt "
                         f"< {spec['min_fs']} pt (use a wider figure or fewer structures per row)")
        x0, y0, x1, y1 = ax.get_position().extents
        if x0 < -1e-3 or y0 < -1e-3 or x1 > 1 + 1e-3 or y1 > 1 + 1e-3:
            probs.append(f"structure {ax.structure_mm[0]:.0f} x {ax.structure_mm[1]:.0f} mm runs off the figure edge")
    if any(getattr(ax, "highlighted", False) for ax in fig.axes):
        notes.append("highlighted groups in structures: say what each colour means in the legend, "
                     "and check the figure still reads in greyscale")
    probs = list(dict.fromkeys(probs))  # one line per distinct problem
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
