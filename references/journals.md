# Journal figure specifications

Values below are what `pubfig.JOURNALS` encodes. Guidelines change — when the user names a specific journal, check its current "Guide for Authors / Artwork" page if you can browse, and override with `pubfig.figsize("generic", width=<mm>)` or `setup(fontsize=...)` if it differs.

| Key | Publisher / examples | Single | 1.5 col | Double | Max height | Font (final size) | Panel label |
|---|---|---|---|---|---|---|---|
| `nature` | Nature, Nat. Commun., Sci. Rep., Nat. family | 89 mm | 120 mm | 183 mm | 170 mm (w/ legend ~247) | 5–7 pt, Arial/Helvetica | **a b c** lowercase bold |
| `science` | Science, Sci. Adv. | 57 mm | 121 mm | 184 mm | ~230 mm | ≥6 pt, Helvetica | **A B C** |
| `cell` | Cell, Cell Rep., Mol. Cell, iScience | 85 mm | 114 mm | 174 mm | ~225 mm | 6–8 pt, Arial | **A B C** |
| `elsevier` | Eur. J. Med. Chem., Bioorg. Chem., J. Mol. Struct., Biomed. Pharmacother., Comput. Biol. Med. | 90 mm | 140 mm | 190 mm | ~240 mm | ≥6 pt (7 pt typical) | (a) or **a** |
| `acs` | J. Med. Chem., J. Chem. Inf. Model., ACS Omega, J. Nat. Prod., Mol. Pharm. | 3.25 in (82.5 mm) | — | 7 in (178 mm) | 9.5 in (240 mm) | ≥4.5 pt (aim 6–8), Arial/Helvetica | **A B C** or (a) |
| `rsc` | RSC Adv., Med. Chem. Commun., PCCP | 83 mm | — | 171 mm | 233 mm | 7–9 pt | (a) |
| `wiley` | ChemMedChem, Chem. Biol. Drug Des., Arch. Pharm. | 80 mm | — | 180 mm | ~230 mm | 7–9 pt | (a) |
| `springer` | J. Comput.-Aided Mol. Des., Mol. Divers., BMC journals | 84 mm | 129 mm | 174 mm | 234 mm | 8–12 pt (Helvetica/Arial) | **a** lowercase bold |
| `plos` | PLOS ONE, PLOS Comput. Biol. | 2.63–7.5 in (5.2 in typical = 132 mm) | — | 7.5 in (190 mm) | 8.75 in (222 mm) | 8–12 pt Arial | **A** |
| `frontiers` | Front. Pharmacol., Front. Chem. | 85 mm | — | 180 mm | — | ≥8 pt | **A** |
| `mdpi` | Molecules, Pharmaceuticals, IJMS | ~85 mm | — | ~170 mm | — | ≥8 pt | (a) |
| `ieee` | IEEE journals | 3.5 in | — | 7.16 in | — | 8–10 pt | (a) |

## Resolution & format

| Content | Minimum DPI at final size | Preferred file |
|---|---|---|
| Line art (plots, charts) | 1000–1200 if raster | Vector: PDF / EPS / SVG |
| Combination (plot + photo/structure image) | 500–600 | TIFF (LZW) or PDF with embedded raster |
| Photographs / micrographs | 300 | TIFF |

- Colour mode: RGB is accepted almost everywhere online; some print journals ask CMYK — convert at the end with the publisher's profile, never design in CMYK (colours shift).
- No alpha channel in TIFF/PNG (`pubfig.save` flattens to white RGB).
- Fonts embedded (TrueType / Type 42). Type 3 fonts trigger failures in Editorial Manager and ScholarOne PDF builders.
- Do not upload JPEG for line art (compression artefacts around text).
- Some journals want each panel as a separate file or a layered source file (AI/EPS) at acceptance; the SVG/PDF output with live text serves that.

## Chemical structures

ACS requires structures in the ChemDraw **ACS Document 1996** style; Elsevier, Wiley, RSC, Springer and MDPI accept it and most medicinal-chemistry papers use it. Some journals (e.g. Angew. Chem., RSC titles) offer their own ChemDraw template: if the user names one, follow it.

| Setting | ACS Document 1996 |
|---|---|
| Fixed bond length | 0.508 cm (14.4 pt) |
| Line width / bold (wedge) width | 0.6 pt / 2.0 pt |
| Double-bond spacing | 18 % of bond length |
| Hash spacing / margin | 2.5 pt / 1.6 pt |
| Atom labels | Arial (or Helvetica) 10 pt |
| Chain angle | 120° |

- Print the structures at 100 % where possible. If the figure must shrink, every structure shrinks by the same factor; keep atom labels ≥ 6 pt (scale ≥ 0.6).
- Same bond length in every structure of a figure, and the common core drawn in the same orientation across a series.
- Kekulé rings (alternating double bonds), not aromatic circles. Show stereochemistry with wedge/hash bonds; define any abbreviations (Boc, Ph, R groups) in the legend.
- Compound numbers in bold, consistent with the text and tables.
- Submit as vector (PDF/EPS) or TIFF at 1200 dpi (line art). Keep the editable source (.cdx / .mol) — journals and reviewers ask for it.

## Image-integrity rules (blots, gels, microscopy)

Q1 journals screen for manipulation. When a figure contains such images:
- Only linear brightness/contrast adjustments, applied to the whole image and to controls equally.
- Cropping must be indicated; keep uncropped originals for Supplementary (many journals require them).
- Splicing lanes requires a visible dividing line and a note in the legend.
- Scale bars on every micrograph (not magnification in text); state what the scale bar represents in the legend.
- Never use the same image to represent different conditions.

## Legends checklist (reviewers look for these)
- Title sentence stating the finding.
- Every panel described; every abbreviation and symbol defined.
- n, replicate type, what centre and error bars are (mean ± SD, median and IQR...).
- Statistical test, correction for multiple comparisons, exact P or star definitions.
- For computational work: software + version, force field, simulation length, number of replicas, docking program + scoring function, structure PDB IDs.
