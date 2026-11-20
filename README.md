# publication-figures

Agent skill for making journal-ready (Q1) scientific figures with Python/matplotlib: exact column widths, embedded TrueType fonts, 600 dpi RGB TIFF, colour-blind-safe palettes, and statistics reporting that reviewers expect.

Plain `SKILL.md` + Python, so it works with any agent that reads skills or instruction files (Claude Code, Codex, Gemini CLI, Cursor, ...).

- **Journals:** Nature, Science, Cell, Elsevier, ACS, RSC, Wiley, Springer, PLOS, Frontiers, MDPI, IEEE, generic
- **Domains:** medicinal chemistry, bioinformatics, molecular dynamics, general
- **Recipes (14):** group comparison + stats, dose-response/IC50 (4PL + CI), MD RMSD/Rg/RMSF/H-bonds (GROMACS `.xvg`), docking, volcano, clustered heatmap, QSAR, Kaplan-Meier, ROC, ADMET, multi-panel layouts, chemical structure series, structure + data SAR panels
- **Chemical structures:** vector drawings in ACS Document 1996 style (bond 14.4 pt, line 0.6 pt, Arial 10 pt) from SMILES/MOL via RDKit, one common scale, scaffold-aligned series, R-group labels, bold compound numbers, `.mol` export for ChemDraw
- **Checks:** `pubfig.audit()` on the live figure (sizes, fonts, structure scale, bond-through-label guard), `pubfig.py check` on saved PDF/TIFF (width, Type 3 fonts, dpi, alpha)

Requires Python 3 with `matplotlib` and `numpy` (`scipy`, `pandas`, `lifelines`, `scikit-learn` only for some recipes; `rdkit` for structures). No install needed: `uv run --with matplotlib --with numpy python script.py`.

## Install

Claude Code (personal skills):

```bash
git clone https://github.com/neuxora/publication-figures ~/.claude/skills/publication-figures
```

Other agents: clone anywhere and point your `AGENTS.md` / `GEMINI.md` / rules file at `SKILL.md`.

## Layout

| Path | Purpose |
|---|---|
| `SKILL.md` | Workflow and Q1 rules the agent follows |
| `scripts/pubfig.py` | Journal presets, helpers, audit and file checker (`python pubfig.py demo OUT`) |
| `references/journals.md` | Per-publisher size/format specs |
| `references/recipes.md` | 12 runnable recipes |
| `evals/` | 3 test cases used to benchmark the skill |

## Benchmark (iteration 1)

3 cases, with vs. without the skill: **26/26** vs. **22/26** checks passed. Baseline failures were all file-spec problems (RGBA TIFFs, a 180 mm instead of 183 mm figure); scientific content was correct in both. Journal specs change; always confirm against the current author guidelines.

## License

MIT
