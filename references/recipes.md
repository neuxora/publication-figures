# Recipes

Each block is self-contained and runs on synthetic data — replace the data section with the user's data. All assume `pubfig.py` is importable and start with:

```python
import numpy as np, matplotlib.pyplot as plt
import pubfig
spec = pubfig.setup("nature")
```

Contents
1. Group comparison + statistics (box/strip, ANOVA/Dunnett or Kruskal–Wallis)
2. Dose–response, 4PL fit, IC50 with 95% CI
3. MD time series (RMSD / Rg / SASA / H-bonds), GROMACS `.xvg` reader
4. RMSF per residue
5. Docking: ranked binding energies + residue-interaction heatmap
6. Volcano plot
7. Clustered heatmap (z-scores)
8. Regression / QSAR predicted vs observed
9. Kaplan–Meier survival
10. ROC curve
11. ADMET / property heatmap (radar alternative)
12. Multi-panel mosaics with mixed sizes
13. Chemical structures: compound series / scaffold with R groups (RDKit, ACS 1996 style)
14. Structures as a panel next to data (SAR figure)

---

## 1. Group comparison + statistics

```python
from scipy import stats
rng = np.random.default_rng(1)
groups = ["Vehicle", "5 µM", "10 µM", "20 µM"]
data = [rng.normal(m, 7, 6) for m in (100, 82, 61, 40)]          # n = 6 biological replicates
colors = [pubfig.GREY] + pubfig.OKABE_ITO[:3]

# check assumptions -> choose test
normal = all(stats.shapiro(d).pvalue > 0.05 for d in data)
equal_var = stats.levene(*data).pvalue > 0.05
if normal and equal_var:
    omni = stats.f_oneway(*data); test = "one-way ANOVA"
    pvals = stats.dunnett(*data[1:], control=data[0]).pvalue     # vs control; stats.tukey_hsd(*data) for all pairs
    post_name = "Dunnett's test vs vehicle"
else:
    omni = stats.kruskal(*data); test = "Kruskal–Wallis"
    raw = [stats.mannwhitneyu(data[0], d).pvalue for d in data[1:]]
    pvals = stats.false_discovery_control(raw, method="bh")
    post_name = "Mann–Whitney U, Benjamini–Hochberg corrected"

fig, ax = plt.subplots(figsize=pubfig.figsize("nature", "single", 0.8))
for i, (d, c) in enumerate(zip(data, colors)):
    ax.scatter(i + rng.uniform(-0.12, 0.12, d.size), d, s=10, color=c, lw=0, zorder=3, alpha=0.9)
    m, sd = d.mean(), d.std(ddof=1)
    ax.errorbar(i, m, yerr=sd, fmt="none", ecolor="k", elinewidth=0.8, capsize=3, zorder=4)
    ax.hlines(m, i - 0.22, i + 0.22, color="k", lw=1.2, zorder=4)
ax.set_xticks(range(len(groups)), groups)
ax.set_xlabel("Treatment")
ax.set_ylabel("Cell viability (% of vehicle)")
ax.set_ylim(0, 150)
ytop = max(d.max() for d in data) * 1.04
for k, p in enumerate(pvals, start=1):                           # brackets stacked upward
    ytop = pubfig.sig_bracket(ax, 0, k, ytop, p=p) + 5
pubfig.audit(fig, "nature", "single")
pubfig.save(fig, "figures/fig_groups")
print(test, omni, post_name, pvals)
```
Legend: "Mean ± SD, n = 6 biological replicates (dots). One-way ANOVA (F = …, P = …) with Dunnett's test vs vehicle."
- Set the y limit after placing brackets so they are not clipped.
- Paired designs (same animal/cell line before–after): `stats.ttest_rel` / `stats.wilcoxon`, and draw lines connecting pairs.
- Two groups only: Welch's t-test (`stats.ttest_ind(a, b, equal_var=False)`) is the safe default.

## 2. Dose–response (4PL) with IC50 and 95% CI

```python
from scipy.optimize import curve_fit
from scipy import stats

def fourpl(logc, bottom, top, logic50, hill):
    return bottom + (top - bottom) / (1 + 10 ** ((logc - logic50) * hill))

conc = np.array([0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100])        # µM
rng = np.random.default_rng(2)
compounds = {"Cpd 1": 0.8, "Cpd 2": 4.0, "Doxorubicin": 0.3}       # true IC50s for the synthetic data
reps = {k: np.array([fourpl(np.log10(conc), 2, 100, np.log10(v), 1.1) + rng.normal(0, 4, conc.size) for _ in range(3)])
        for k, v in compounds.items()}                            # shape (3 replicates, n conc)

fig, ax = plt.subplots(figsize=pubfig.figsize("acs", "single", 0.8))
xfit = np.logspace(np.log10(conc.min()) - 0.3, np.log10(conc.max()) + 0.3, 300)
colors = dict(zip(compounds, pubfig.OKABE_ITO[:2] + ["#000000"]))  # reference drug in black
results = {}
for name, y in reps.items():
    X = np.tile(np.log10(conc), y.shape[0]); Y = y.ravel()       # fit ALL replicates, not means
    p0 = [Y.min(), Y.max(), np.median(np.log10(conc)), 1.0]
    popt, pcov = curve_fit(fourpl, X, Y, p0=p0, maxfev=20000)
    se = np.sqrt(np.diag(pcov))[2]
    tcrit = stats.t.ppf(0.975, len(Y) - 4)
    ic50, lo, hi = 10 ** popt[2], 10 ** (popt[2] - tcrit * se), 10 ** (popt[2] + tcrit * se)
    results[name] = (ic50, lo, hi, popt[3])
    c = colors[name]
    ax.errorbar(conc, y.mean(0), yerr=y.std(0, ddof=1), fmt="o", ms=3.5, color=c, mec="none",
                elinewidth=0.7, capsize=2)
    ax.plot(xfit, fourpl(np.log10(xfit), *popt), color=c, lw=1,
            label=f"{name} (IC$_{{50}}$ = {ic50:.2f} µM)")
ax.set_xscale("log")
ax.set_xlabel("Concentration (µM)")
ax.set_ylabel("Cell viability (%)")
ax.set_ylim(-5, 120)
ax.legend(loc="lower left")
pubfig.audit(fig, "acs", "single")
pubfig.save(fig, "figures/fig_dose_response")
for k, (ic, lo, hi, h) in results.items():
    print(f"{k}: IC50 = {ic:.2f} µM (95% CI {lo:.2f}–{hi:.2f}), Hill = {h:.2f}")
```
- Report IC50 with 95% CI (and pIC50 ± SE if the field prefers) in a table or the legend.
- If the curve doesn't reach the bottom plateau, say "IC50 > highest tested concentration" rather than extrapolating.
- Normalise to vehicle (100%) and blank/positive control (0%) before fitting if raw absorbances are given.

## 3. MD time series (RMSD, Rg, SASA, H-bonds)

```python
def read_xvg(path):
    """GROMACS .xvg -> (array[n, k], legends). Skips # and @ lines; time column usually in ps."""
    rows, legends = [], []
    for line in open(path):
        if line.startswith("@") and " legend " in line and line.split()[1].startswith("s"):
            legends.append(line.split('"')[1])
        elif line.strip() and line[0] not in "#@":
            rows.append([float(v) for v in line.split()])
    return np.array(rows), legends

def running_mean(y, w):
    return np.convolve(y, np.ones(w) / w, mode="valid")

# synthetic: 3 systems, 100 ns, 1 frame / 10 ps
t = np.arange(0, 100.01, 0.01)                                      # ns  (xvg time in ps -> /1000)
rng = np.random.default_rng(3)
systems = {"Apo": 0.30, "Complex 1": 0.20, "Complex 2": 0.25}
rmsd = {k: v * (1 - np.exp(-t / 5)) + 0.02 * rng.standard_normal(t.size) + 0.05 for k, v in systems.items()}
rg = {k: 1.85 + 0.02 * rng.standard_normal(t.size) + 0.01 * i for i, k in enumerate(systems)}
colors = dict(zip(systems, [pubfig.GREY] + pubfig.OKABE_ITO[:2]))

fig, axs = plt.subplots(1, 2, figsize=pubfig.figsize("elsevier", "double", 0.36))
w = 100                                                              # 1 ns window
for name in systems:
    for ax, series in ((axs[0], rmsd), (axs[1], rg)):
        ax.plot(t, series[name], color=colors[name], lw=0.4, alpha=0.35)
        ax.plot(t[w - 1:], running_mean(series[name], w), color=colors[name], lw=1, label=name)
axs[0].set_ylabel("Backbone RMSD (nm)")
axs[1].set_ylabel("Radius of gyration (nm)")
for ax in axs:
    ax.set_xlabel("Time (ns)")
    ax.set_xlim(0, 100)
axs[0].legend(loc="lower right")
pubfig.label_panels(axs, pubfig.spec_for("elsevier"))
pubfig.audit(fig, "elsevier", "double")
pubfig.save(fig, "figures/fig_md")
for k in systems:   # report equilibrated averages (e.g. last 50 ns) in text/table
    m = t >= 50; print(k, f"RMSD {rmsd[k][m].mean():.3f} ± {rmsd[k][m].std():.3f} nm")
```
- Replicate simulations: plot the mean line with a shaded ±SD band (`ax.fill_between`) instead of overlaying all replicas.
- Add a marginal histogram/KDE (narrow right-hand axes, `sharey`) when the distribution matters.
- Units: GROMACS outputs nm and ps; AMBER/Desmond often Å — label the actual unit. Convert ps → ns.
- H-bond counts are integers: plot a running mean or a histogram, not a raw jagged line.
- Raw lines with thousands of points: `rasterized=True` keeps the PDF light.

## 4. RMSF per residue

```python
rng = np.random.default_rng(4)
res = np.arange(1, 301)
base = 0.08 + 0.05 * np.abs(np.sin(res / 15)) + np.where((res > 140) & (res < 160), 0.25, 0)
rmsf = {"Apo": base + 0.02 * rng.random(res.size), "Complex 1": base * 0.8 + 0.02 * rng.random(res.size)}
fig, ax = plt.subplots(figsize=pubfig.figsize("nature", "double", 0.28))
for (k, v), c in zip(rmsf.items(), [pubfig.GREY, pubfig.OKABE_ITO[0]]):
    ax.plot(res, v, lw=0.8, color=c, label=k)
ax.axvspan(140, 160, color="0.9", lw=0, zorder=0)                    # e.g. activation loop
ax.text(150, 0.42, "A-loop", ha="center", va="bottom")
for r in (45, 98, 210):                                             # binding-site residues
    ax.axvline(r, color=pubfig.OKABE_ITO[3], lw=0.5, ls="--")
ax.set_xlabel("Residue number"); ax.set_ylabel("RMSF (nm)"); ax.set_xlim(res.min(), res.max())
ax.legend(loc="upper right")
pubfig.audit(fig, "nature", "double")
pubfig.save(fig, "figures/fig_rmsf")
```

## 5. Docking: ranking + interaction heatmap

```python
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Patch
ligs = np.array(["Ref (co-crystal)", "L1", "L2", "L3", "L4", "L5", "L6"])
dg = np.array([-8.4, -9.6, -9.1, -8.8, -7.9, -7.2, -6.5])          # kcal/mol (more negative = better)
order = np.argsort(dg)[::-1]                                        # barh draws bottom-up -> best at top
residues = ["Asp86", "Lys33", "Leu83", "Glu81", "Phe80", "Ile10"]
inter = np.random.default_rng(5).integers(0, 4, (len(ligs), len(residues)))  # 0 none, 1 hydrophobic, 2 H-bond, 3 both

fig, axs = plt.subplots(1, 2, figsize=pubfig.figsize("acs", "double", 0.38), width_ratios=[1, 1.2])
colors = ["#000000" if l.startswith("Ref") else pubfig.OKABE_ITO[0] for l in ligs[order]]
axs[0].barh(ligs[order], dg[order], color=colors, height=0.6)
axs[0].axvline(dg[0], color="k", lw=0.6, ls="--")                   # reference ligand score
axs[0].set_xlabel("Binding affinity (kcal mol$^{-1}$)")
axs[0].set_ylabel("Ligand")
cmap = ListedColormap(["#FFFFFF", pubfig.OKABE_ITO[1], pubfig.OKABE_ITO[0], pubfig.OKABE_ITO[3]])
axs[1].imshow(inter, cmap=cmap, norm=BoundaryNorm(np.arange(-0.5, 4), 4), aspect="auto")
axs[1].set_xticks(range(len(residues)), residues, rotation=45, ha="right")
axs[1].set_yticks(range(len(ligs)), ligs)
axs[1].set_xticks(np.arange(-.5, len(residues)), minor=True); axs[1].set_yticks(np.arange(-.5, len(ligs)), minor=True)
axs[1].grid(which="minor", color="0.8", lw=0.4); axs[1].tick_params(which="minor", length=0)
axs[1].set_xlabel("Residue")
axs[1].legend([Patch(fc=c, ec="0.6", lw=0.4) for c in cmap.colors[1:]], ["Hydrophobic", "H-bond", "Both"],
              loc="upper center", bbox_to_anchor=(0.5, -0.35), ncol=3)
pubfig.label_panels(axs, pubfig.spec_for("acs"))
pubfig.audit(fig, "acs", "double")
pubfig.save(fig, "figures/fig_docking")
```
- State docking program, scoring function, grid box, PDB ID and redocking RMSD of the co-crystal ligand (validation) in the legend/methods.
- Docking scores differing by < ~1 kcal/mol are within scoring error; don't add significance stars to docking scores.

## 6. Volcano plot

```python
rng = np.random.default_rng(6)
n = 5000
lfc = rng.normal(0, 1.2, n); padj = np.clip(10 ** (-np.abs(lfc) * rng.uniform(0.5, 3, n)), 1e-50, 1)
genes = np.array([f"G{i}" for i in range(n)])
fc_t, p_t = 1.0, 0.05
up = (lfc >= fc_t) & (padj < p_t); down = (lfc <= -fc_t) & (padj < p_t); ns = ~(up | down)
y = -np.log10(padj)

fig, ax = plt.subplots(figsize=pubfig.figsize("nature", "single", 0.9))
ax.scatter(lfc[ns], y[ns], s=3, c="0.75", lw=0, rasterized=True, label=f"NS ({ns.sum()})")
ax.scatter(lfc[up], y[up], s=4, c=pubfig.UP, lw=0, rasterized=True, label=f"Up ({up.sum()})")
ax.scatter(lfc[down], y[down], s=4, c=pubfig.DOWN, lw=0, rasterized=True, label=f"Down ({down.sum()})")
for v in (-fc_t, fc_t): ax.axvline(v, color="0.4", lw=0.5, ls="--")
ax.axhline(-np.log10(p_t), color="0.4", lw=0.5, ls="--")
top = np.argsort(np.where(ns, np.inf, padj))[:8]                    # label top significant hits
for i in top:
    ax.annotate(genes[i], (lfc[i], y[i]), xytext=(3, 2), textcoords="offset points", fontsize=spec["fs"] - 1)
ax.set_xlabel("log$_2$ fold change"); ax.set_ylabel("−log$_{10}$ adjusted P")
ax.legend(loc="upper left", markerscale=2, handletextpad=0.2)
pubfig.audit(fig, "nature", "single")
pubfig.save(fig, "figures/fig_volcano")
```
- `rasterized=True` on thousands of points keeps the PDF small while text/axes stay vector.
- Label collisions: `pip install adjustText` then `adjust_text(texts, ax=ax)`; otherwise label fewer genes.
- Always adjusted P (BH/FDR) on the y axis for omics; state thresholds in the legend.

## 7. Clustered heatmap (z-scores)

```python
from scipy.cluster.hierarchy import linkage, leaves_list, dendrogram
from matplotlib.colors import TwoSlopeNorm
rng = np.random.default_rng(7)
mat = rng.normal(0, 1, (30, 8)); mat[:10, 4:] += 2; mat[10:20, :4] += 2
rows = np.array([f"Gene{i}" for i in range(30)]); cols = np.array([f"S{j}" for j in range(8)])
z = (mat - mat.mean(1, keepdims=True)) / mat.std(1, ddof=1, keepdims=True)   # row z-score
rl = linkage(z, "average", metric="correlation"); ro = leaves_list(rl)
cl = linkage(z.T, "average", metric="correlation"); co = leaves_list(cl)

fig = plt.figure(figsize=pubfig.figsize("nature", "single", 1.3))
gs = fig.add_gridspec(2, 3, width_ratios=[0.18, 1, 0.05], height_ratios=[0.12, 1])
axd = fig.add_subplot(gs[1, 0]); axt = fig.add_subplot(gs[0, 1]); ax = fig.add_subplot(gs[1, 1]); cax = fig.add_subplot(gs[1, 2])
dendrogram(rl, ax=axd, orientation="left", no_labels=True, color_threshold=0, above_threshold_color="k")
dendrogram(cl, ax=axt, no_labels=True, color_threshold=0, above_threshold_color="k")
for a in (axd, axt): a.axis("off")
axd.invert_yaxis()                                                  # match imshow row order
lim = np.abs(z).max()
im = ax.imshow(z[ro][:, co], cmap="RdBu_r", norm=TwoSlopeNorm(0, -lim, lim), aspect="auto", interpolation="none")
ax.set_xticks(range(len(cols)), cols[co], rotation=90)
ax.set_yticks(range(len(rows)), rows[ro]); ax.yaxis.tick_right()
ax.tick_params(length=0)
for s in ax.spines.values(): s.set_visible(False)
cb = fig.colorbar(im, cax=cax); cb.set_label("Row z-score"); cb.outline.set_linewidth(0.5)
ax.set_xlabel("Sample")
pubfig.audit(fig, "nature", "single")
pubfig.save(fig, "figures/fig_heatmap")
```
- With many rows, drop row labels and label only clusters; tick labels must stay ≥ journal minimum.
- Correlation matrices: `RdBu_r` with fixed limits −1..1, lower triangle only, annotate values if ≤ ~10×10.

## 8. Regression / QSAR predicted vs observed

```python
rng = np.random.default_rng(8)
obs_tr = rng.uniform(4, 9, 80); pred_tr = obs_tr + rng.normal(0, 0.4, 80)
obs_te = rng.uniform(4, 9, 20); pred_te = obs_te + rng.normal(0, 0.55, 20)
def metrics(o, p):
    r2 = 1 - np.sum((o - p) ** 2) / np.sum((o - o.mean()) ** 2); rmse = np.sqrt(np.mean((o - p) ** 2))
    return r2, rmse

fig, ax = plt.subplots(figsize=pubfig.figsize("acs", "single", 0.95))
lims = [3.5, 9.5]
ax.plot(lims, lims, color="k", lw=0.7, ls="--", zorder=1)           # y = x
for o, p, c, mk, lab in ((obs_tr, pred_tr, pubfig.OKABE_ITO[0], "o", "Training"),
                         (obs_te, pred_te, pubfig.OKABE_ITO[3], "^", "Test")):
    r2, rmse = metrics(o, p)
    ax.scatter(o, p, s=14, color=c, marker=mk, lw=0.3, ec="white", zorder=3,
               label=f"{lab} (n = {o.size}, R$^2$ = {r2:.2f}, RMSE = {rmse:.2f})")
ax.set_xlim(lims); ax.set_ylim(lims); ax.set_aspect("equal")
ax.set_xlabel("Experimental pIC$_{50}$"); ax.set_ylabel("Predicted pIC$_{50}$")
ax.legend(loc="upper left")
pubfig.audit(fig, "acs", "single")
pubfig.save(fig, "figures/fig_qsar")
```
Simple correlation instead: `res = scipy.stats.linregress(x, y)`; draw the line plus a 95% CI band of the fit; put `r = …, P = …, n = …` in a corner.

## 9. Kaplan–Meier

```python
# pip install lifelines
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
rng = np.random.default_rng(9)
tA, tB = rng.exponential(20, 40), rng.exponential(35, 40); eA, eB = rng.random(40) < 0.8, rng.random(40) < 0.7
fig, ax = plt.subplots(figsize=pubfig.figsize("nature", "single", 0.8))
for t, e, lab, c in ((tA, eA, "Vehicle", pubfig.GREY), (tB, eB, "Treated", pubfig.OKABE_ITO[0])):
    KaplanMeierFitter().fit(t, e, label=lab).plot_survival_function(
        ax=ax, ci_show=True, color=c, lw=1, show_censors=True, censor_styles={"ms": 4, "marker": "|"})
p = logrank_test(tA, tB, eA, eB).p_value
ax.text(0.97, 0.95, f"Log-rank {pubfig.format_p(p)}", transform=ax.transAxes, ha="right", va="top")
ax.set_xlabel("Time (days)"); ax.set_ylabel("Survival probability"); ax.set_ylim(0, 1.02)
pubfig.audit(fig, "nature", "single")
pubfig.save(fig, "figures/fig_km")
```
Clinical journals expect a numbers-at-risk table under the axis (`lifelines.plotting.add_at_risk_counts`).

## 10. ROC

```python
# pip install scikit-learn
from sklearn.metrics import roc_curve, roc_auc_score
rng = np.random.default_rng(10)
y = rng.integers(0, 2, 300); s1 = y + rng.normal(0, 0.7, 300); s2 = y + rng.normal(0, 1.2, 300)
fig, ax = plt.subplots(figsize=pubfig.figsize("nature", "single", 1.0))
ax.plot([0, 1], [0, 1], color="0.6", lw=0.6, ls="--")
for s, lab, c in ((s1, "Model A", pubfig.OKABE_ITO[0]), (s2, "Model B", pubfig.OKABE_ITO[1])):
    fpr, tpr, _ = roc_curve(y, s)
    ax.plot(fpr, tpr, color=c, lw=1, label=f"{lab} (AUC = {roc_auc_score(y, s):.2f})")
ax.set_aspect("equal"); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
ax.set_xlabel("False positive rate"); ax.set_ylabel("True positive rate"); ax.legend(loc="lower right")
pubfig.audit(fig, "nature", "single")
pubfig.save(fig, "figures/fig_roc")
```
Report AUC with a 95% CI (bootstrap ~1000×) when comparing models.

## 11. ADMET / property heatmap (radar alternative)

```python
props = ["MW", "cLogP", "TPSA", "HBD", "HBA", "RotB"]
limits = np.array([500, 5, 140, 5, 10, 10])                        # Lipinski / Veber limits
cpds = ["Cpd 1", "Cpd 2", "Cpd 3", "Ref"]
vals = np.array([[412, 3.1, 88, 2, 6, 5], [530, 5.6, 72, 1, 7, 8], [355, 2.2, 120, 3, 8, 4], [394, 4.0, 95, 2, 6, 6]])
frac = vals / limits                                                # 1.0 = at the limit
fig, ax = plt.subplots(figsize=pubfig.figsize("acs", "single", 0.6))
im = ax.imshow(frac, cmap="cividis_r", vmin=0, vmax=1.3, aspect="auto")
for i in range(len(cpds)):
    for j in range(len(props)):
        viol = frac[i, j] > 1
        ax.text(j, i, f"{vals[i, j]:g}" + ("*" if viol else ""), ha="center", va="center",
                color="white" if frac[i, j] > 0.75 else "black", fontweight="bold" if viol else "normal")
ax.set_xticks(range(len(props)), props); ax.set_yticks(range(len(cpds)), cpds)
ax.set_xlabel("Property"); ax.tick_params(length=0)
for s in ax.spines.values(): s.set_visible(False)
cb = fig.colorbar(im, ax=ax, shrink=0.9); cb.set_label("Fraction of rule limit")
pubfig.audit(fig, "acs", "single")
pubfig.save(fig, "figures/fig_admet")
```
If the user insists on a radar (SwissADME-style bioavailability radar): `subplot_kw={"projection": "polar"}`, normalise each axis to its limit, draw the acceptable range as a grey band, ≤ 4 compounds per radar.

## 12. Multi-panel mosaics

```python
fig, axd = plt.subplot_mosaic([["a", "a", "b"], ["c", "d", "e"]],
                              figsize=pubfig.figsize("nature", "double", 0.55),
                              height_ratios=[1.2, 1])
for k, ax in axd.items():                                          # ... real plotting here
    ax.plot([0, 1], [0, 1]); ax.set_xlabel("x (unit)"); ax.set_ylabel("y (unit)")
pubfig.label_panels([axd[k] for k in "abcde"], spec)
pubfig.audit(fig, "nature", "double")
pubfig.save(fig, "figures/fig_mosaic")
```
- `sharex/sharey` (or `ax.sharey(other)`) for panels compared side by side.
- Chemical structures: draw them as vectors with `pubfig.draw_structures` (recipes 13–14). Image panel (docking pose from PyMOL, microscopy): `ax.imshow(plt.imread("pose.png")); ax.axis("off")` — export that PNG at ≥ 600 dpi at its printed size.
- One shared legend for all panels: `fig.legend(handles, labels, loc="outside upper center", ncol=n)` (works with constrained layout).

## 13. Chemical structures (compound series, scaffold with R groups)

Needs `pip install rdkit`. Structures are drawn as vectors in **ACS Document 1996** style — the ChemDraw preset that ACS requires and Elsevier, Wiley, RSC, MDPI and Springer accept: bond 14.4 pt (0.508 cm), line 0.6 pt, bold/wedge 2.0 pt, double-bond spacing 18 %, Arial 10 pt atom labels, Kekulé rings, wedge/hash stereo.

```python
from rdkit import Chem
spec = pubfig.setup("acs")
series = {"4a": "O=C(Nc1ccc(Cl)cc1)c1ccc2[nH]ccc2c1",          # from the user's SMILES / SDF
          "4b": "O=C(Nc1ccc(OC)cc1)c1ccc2[nH]ccc2c1",
          "4c": "O=C(Nc1ccc(C(F)(F)F)cc1)c1ccc2[nH]ccc2c1"}
scaffold = "O=C(Nc1ccc([*:1])cc1)c1ccc2n([*:2])ccc2c1"         # [*:n] is drawn as R with subscript n

fig = plt.figure(figsize=pubfig.figsize("acs", "double", height_mm=45))
axs, scale = pubfig.draw_structures(fig, [scaffold, *series.values()], ["General structure", *series],
                                    core="O=C(N)c1ccc2[nH]ccc2c1")   # same scaffold orientation in every compound
print(f"structure scale {scale:.2f} (1.00 = ACS 1996 at 100 %)")
pubfig.audit(fig, "acs", "double")
pubfig.save(fig, "figures/fig_structures", dpi=1200)               # line art: 1200 dpi TIFF
for lab, smi in series.items():                                    # editable source: opens in ChemDraw/Marvin
    Chem.MolToMolFile(Chem.MolFromSmiles(smi), f"figures/{lab}.mol")
```
- `draw_structures` uses **one scale for all structures** (equal bond lengths) and the largest that fits, never above 100 %. If the scale drops below ~0.6, atom labels go under 6 pt and `audit` flags it: use fewer per row (`ncols=`), a taller or wider figure, or split into two figures.
- `core=` aligns the shared scaffold (SMILES or SMARTS) so a series reads at a glance; compounds without it keep their own layout.
- Compound numbers are bold (ACS/RSC/Elsevier convention); keep the same numbers as the text.
- Already have ChemDraw/Marvin drawings? Load the MOL file (`Chem.MolFromMolFile("4a.mol")`) — its 2D coordinates are kept, so the author's orientation survives.
- `mol_svg` warns if RDKit draws a bond through an atom label; look at the render, or give it a MOL file with hand-placed coordinates.
- Black and white only. Highlight a change with a coloured atom/bond only if the legend explains it.

## 14. Structures as a panel next to data (SAR figure)

```python
spec = pubfig.setup("acs")
names = ["4a", "4b", "4c"]
smiles = ["O=C(Nc1ccc(Cl)cc1)c1ccc2[nH]ccc2c1", "O=C(Nc1ccc(OC)cc1)c1ccc2[nH]ccc2c1",
          "O=C(Nc1ccc(C(F)(F)F)cc1)c1ccc2[nH]ccc2c1"]
pic50, lo, hi = np.array([5.50, 4.81, 4.0]), np.array([5.45, 4.71, np.nan]), np.array([5.55, 4.91, np.nan])

fig, ax = plt.subplot_mosaic([["s", "s", "b"]], figsize=pubfig.figsize("acs", "double", height_mm=60))
b = ax["b"]
for i, (p, l, h) in enumerate(zip(pic50, lo, hi)):
    if np.isnan(l):                                                   # inactive: IC50 > 100 µM
        b.plot(i, p, marker="v", color=pubfig.GREY, ms=5, ls="none")
    else:
        b.errorbar(i, p, yerr=[[p - l], [h - p]], fmt="o", color=pubfig.OKABE_ITO[i], ms=4, capsize=2, lw=0.8)
b.set_xticks(range(3), names); [t.set_fontweight("bold") for t in b.get_xticklabels()]
b.set_xlim(-0.6, 2.6); b.set_ylabel("pIC$_{50}$ (MCF-7)")
pubfig.draw_structures(fig, smiles, names, rect_mm=ax["s"], core="O=C(N)c1ccc2[nH]ccc2c1")  # fills panel a
pubfig.label_panels([ax["s"], ax["b"]], spec)
pubfig.audit(fig, "acs", "double")
pubfig.save(fig, "figures/fig_sar", dpi=1200)
```
- Pass an empty panel axes as `rect_mm`: the structures fill it after constrained layout has placed the other panels.
- Use the same compound order and numbers in both panels; legend: "Inactive compounds (IC50 > 100 µM) are shown at the assay limit (▼)".
- For interaction diagrams (ligand–residue contacts) use the docking program's 2D diagram (Discovery Studio, LigPlot+, PoseView) exported as vector or ≥ 600 dpi, placed as an image panel (recipe 12).
