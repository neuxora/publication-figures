---
name: publication-figures
description: Make journal-ready scientific figures (Q1 / high-impact quality) with Python matplotlib - exact column widths, correct fonts and DPI, colour-blind-safe palettes, proper statistics annotation, multi-panel layout, vector + TIFF export, and a figure legend draft. Use this skill whenever the user wants a chart, graph, plot, "grafik", "gambar", or "figure" for a paper, thesis, manuscript, journal submission, poster, or reviewer revision - including dose-response/IC50, docking scores, molecular dynamics (RMSD, RMSF, Rg, SASA, H-bonds), ADMET, QSAR, box/violin plots with significance tests, volcano plots, heatmaps, survival curves, ROC, regression, chemical structures ("struktur kimia": compound series, scaffolds with R groups, SAR panels in ACS 1996 style), or any multi-panel figure - and also when they ask to fix, restyle, or "make publication quality" an existing plot, even if they never say "publication".
---

# Publication figures (Q1 standard)

Goal: a figure an editor at Nature / Cell / ACS / Elsevier would accept without a "please revise figures" note, plus the code that reproduces it. Reviewers judge figures fast; most rejection-level problems are mechanical (unreadable 4-pt text, rainbow colours, bars hiding n=3, undefined error bars, blurry raster). This skill makes those mechanical things automatic so effort goes into choosing the right plot.

Works in any agent that can run Python. Needs `matplotlib` + `numpy`; `scipy`/`pandas` only for fitting/stats recipes, `rdkit` only for chemical structures. If a library is missing, install it (`pip install matplotlib scipy pandas`) or run with `uv run --with matplotlib --with scipy --with pandas python script.py`.

## Workflow

### 1. Pin down the brief (ask only what you can't infer)
- **Target journal / publisher** → sets width, font size, panel-label style. Unknown → `generic` and say so.
- **Width**: single column, 1.5, or double. Multi-panel (≥3 panels) usually = double.
- **The message**: one sentence the figure must prove ("Compound 3 is 10× more potent than the reference"). Every design choice serves this sentence.
- **Data shape**: n per group, replicate type (technical vs biological), what the error bars should be.

If the user gave data inline or as a file, read it fully and check it (units, missing values, obvious outliers) before plotting. Don't invent data; if something needed is missing (e.g. no replicate values, only means), say so and plot what is defensible.

### 2. Choose the plot for the data, not out of habit
| Data | Use | Avoid |
|---|---|---|
| Groups, small n (<~15) | Points (jittered) + mean±SD or median line; box+points | Bar+SEM alone (hides n and distribution) |
| Groups, larger n | Box/violin + points (subsample if huge) | Dynamite bars |
| Dose / concentration–response | Points mean±SD on log x + fitted 4PL curve, IC50/EC50 with 95% CI | Linear x axis, connecting dots |
| Time series (MD, kinetics) | Thin raw line (α≈0.3–0.5) + running average; shaded SD for replicates | Thick jagged lines, markers on 5000 points |
| Per-residue (RMSF, contacts) | Line/step vs residue number, annotate key residues/regions | Bars per residue |
| Ranking (docking scores, importance) | Horizontal sorted dot/bar plot, reference highlighted | Unsorted vertical bars, 3D |
| Two continuous variables | Scatter + fit + 95% CI band + r/R², p, n; y=x line for predicted-vs-observed | Unlabelled fit |
| Matrix (expression, interactions) | Heatmap, perceptually uniform/diverging cmap, clustered if meaningful | jet/rainbow |
| Omics differential | Volcano (log2FC vs −log10 adj. P) with thresholds | Unadjusted P |
| Proportions | Stacked bar or dot plot | Pie / donut, 3D anything |
| Many properties per compound (ADMET) | Heatmap of normalized properties, or small multiples | Radar unless the user insists (area distorts) |
| Chemical structures (series, scaffold + R groups, SAR) | `pubfig.draw_structures`: vector, ACS 1996 style, one scale, scaffold aligned, bold compound numbers | Screenshots/JPEG of ChemDraw, mixed bond lengths, aromatic circles |

Recipes with working code for all of the above are in `references/recipes.md` — read the section you need instead of writing from scratch.

### 3. Build it with `scripts/pubfig.py`
Copy `scripts/pubfig.py` next to the figure script (so the user's script runs standalone and can go in the Supplementary / code repository), then:

```python
import matplotlib.pyplot as plt, numpy as np
import pubfig

spec = pubfig.setup("acs")    # nature, science, cell, elsevier, acs, rsc, wiley, springer, plos, frontiers, mdpi, ieee, generic
fig, axs = plt.subplots(1, 2, figsize=pubfig.figsize("acs", "double", aspect=0.42))
# ... plot, using pubfig.OKABE_ITO colours (already the default cycle), pubfig.GREY for controls
pubfig.label_panels(axs, spec)           # a/b or A/B per journal style
pubfig.audit(fig, "acs", "double")       # fix everything it reports, re-run until clean
pubfig.save(fig, "figures/fig2_dose_response")   # .pdf .svg .tiff(600 dpi, RGB) .png
```

What `setup()` gives you, so you don't fight it: Arial/Helvetica sans at journal size, top/right spines off, outward ticks, frameless legend, Okabe–Ito colour cycle, viridis default cmap, TrueType font embedding (no Type 3 — many submission systems reject those), SVG text kept editable, constrained layout on. Do not call `tight_layout()` or save with `bbox_inches="tight"`: both change the final size and break the exact column width.

Other helpers: `sig_bracket(ax, x1, x2, y, p=..., exact=False)`, `p_to_stars(p)`, `format_p(p)`, `check_file(path, journal, width)`. Structures: `draw_structures(fig, smiles_or_mols, labels, ncols=, rect_mm=<area or empty panel axes>, core=<scaffold>)`, `draw_structure(...)` for one, `mol_svg(mol)` for a standalone ACS-1996 SVG. Run `python pubfig.py check fig.pdf fig.tiff --journal acs --width double` on outputs as the final gate.

### 4. Rules that separate Q1 from "revise figures"
These matter because the figure is printed at its final size, read by people with colour-vision deficiency (~8% of men), and scrutinised for statistical honesty.

**Size and type**
- Design at final print size (that's what `figsize` does). Text 5–8 pt at final size, never below the journal minimum; axis labels ≥ tick labels; one font family throughout.
- Axis labels: quantity + unit in parentheses — `Concentration (µM)`, `RMSD (nm)`, `ΔG (kcal mol⁻¹)`. Use real symbols (µ, °, ⁻¹, −) via unicode or mathtext, not "uM" or "-1".
- Mathtext sub/superscripts (`IC$_{50}$`, `mL$^{-1}$`) print at 70 % of the label size; `audit` prints a `note:` when that falls below the journal minimum. Most journals accept it, but for a strict minimum (ACS 4.5 pt) give those labels the size the note suggests.
- Lines ≥ 0.5 pt, markers big enough to see when printed (~3–4 pt).

**Colour**
- Okabe–Ito for categories (default cycle); grey for control/vehicle/non-significant; the same compound/condition gets the same colour in every panel and every figure of the paper.
- Sequential data → `viridis`/`cividis`; data centred on zero (fold change, z-score, correlation) → diverging `RdBu_r` with `TwoSlopeNorm(vcenter=0)` or symmetric limits. Never `jet`/`rainbow`.
- Don't rely on colour alone when it's cheap to add a second cue (marker shape, line style, direct labels). Red–green pairs are out.

**Statistics (state it or reviewers will ask)**
- Show individual data points whenever n is small. Say what the error bars are (SD, SEM, 95% CI) — SD for describing variability, CI for inference; SEM alone is discouraged by many journals.
- Report n and whether replicates are biological or technical.
- Name the test, check its assumptions (normality, equal variance, paired vs independent), correct for multiple comparisons (Tukey/Dunnett/Holm/BH). Prefer exact P values in the legend or on the plot when the journal does (`format_p`); stars are fine if defined.
- Never fit or test on means only when raw replicates exist. Don't truncate a bar axis above zero.

**Chemical structures** (see `references/journals.md` → Chemical structures)
- ACS Document 1996 style at 100 % (bond 14.4 pt, line 0.6 pt, Arial 10 pt) — `draw_structures` does this; never paste screenshots or JPEGs.
- One scale for every structure in the figure (equal bond lengths); the shared scaffold in the same orientation (`core=`); bold compound numbers that match the text; black and white.
- Structures from the user's SMILES/SDF/MOL only: never draw a structure from a name you "remember" without telling the user to verify it. Prefer their MOL files (keeps their orientation); SMILES otherwise.
- Colour only to make a point (varied substituent, pharmacophore, bioisostere): `highlight={SMARTS: colour}`, `highlight_style="halo"` (default, greyscale-safe) or `"color"`; ≤ 3 Okabe–Ito colours, one meaning each across the paper, stated in the legend. No CPK colouring of every heteroatom in scheme/table structures.
- Export `.mol` files alongside so the user can edit in ChemDraw; save structure figures at 1200 dpi (`pubfig.save(..., dpi=1200)`).

**Layout**
- Panels share axes/limits when they should be compared; align them; label panels with bold letters in the journal's style (Nature/Elsevier lowercase **a**, Cell/ACS/Science uppercase **A**).
- Legends inside the axes' free area or direct labels at line ends; no titles inside the plot (the legend/caption carries the title); no gridlines unless reading exact values matters (then light grey, 0.3 pt).
- Remove chartjunk: no 3D, shadows, gradient fills, boxed legends, redundant colour bars.

### 5. Verify before you hand it over
1. `pubfig.audit(...)` returns no issues.
2. `python pubfig.py check <files> --journal J --width W` → OK for the PDF and the TIFF.
3. **Look at the PNG** (open/view the image). Check for: overlapping labels/legend on data, clipped text, panel letters colliding with tick labels, unreadable dense points, colours consistent across panels. Fix and re-render — the first render is rarely final.

### 6. Deliver
Put outputs in a `figures/` folder (or where the user asked):
- `figX_name.pdf` (+ `.svg`) — vector, for submission/editing in Illustrator/Inkscape
- `figX_name.tiff` — 600 dpi RGB (300 for photos, 1000–1200 for pure line art if the journal requires)
- `figX_name.png` — preview
- `figX_name.py` + `pubfig.py` — reproducible script
- `figX_name_legend.md` — draft figure legend, in the user's language unless the journal language differs:

```
**Figure X. <One-sentence title stating the finding.>** (a) <What is plotted, axes, conditions.> (b) ...
Data are mean ± SD of n = 3 independent experiments (biological replicates), each in technical triplicate; dots show individual replicates.
Curves are four-parameter logistic fits; IC50 values (95% CI) are given in the legend. Statistical significance: one-way ANOVA with Dunnett's post hoc test versus vehicle; *P < 0.05, **P < 0.01, ***P < 0.001.
```

Then tell the user briefly: what each panel shows, the key numbers you computed (IC50s, P values, R²), any judgement call you made (e.g. "used SD not SEM", "log-transformed concentration"), and anything they must confirm (journal width, replicate type).

## Reference files
- `references/recipes.md` — tested code for each plot type: group comparison + stats, dose–response 4PL/IC50, MD time series (incl. GROMACS `.xvg` reader), RMSF, docking ranking & interaction heatmap, volcano, clustered heatmap, regression / QSAR predicted-vs-observed, Kaplan–Meier, ROC, ADMET heatmap, multi-panel mosaics, chemical structure series / scaffolds with R groups, structure + data SAR panels.
- `references/journals.md` — width/height/font/DPI/format table per publisher, and the submission-format gotchas (CMYK, layered files, image-integrity rules for blots/micrographs).
