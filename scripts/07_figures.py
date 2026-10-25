"""All figures in figures/, drawn from committed results and the (uncommitted) work matrices."""
from __future__ import annotations

import sys
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from lifelines import KaplanMeierFitter  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from cos_sig.orthologs import map_signature  # noqa: E402
from cos_sig.paths import FIGURES, RESULTS, data_dir, work_dir  # noqa: E402
from cos_sig.scoring import log_zscore, signature_genes_present  # noqa: E402
from cos_sig.signatures import PUBLISHED_ISA, published_gs  # noqa: E402
from cos_sig.survival import align_survival, logrank_p  # noqa: E402

warnings.filterwarnings("ignore")
D, W = data_dir(), work_dir()
FIGURES.mkdir(exist_ok=True)
INK, INK2, MUTED, SURF = "#0b0b0b", "#52514e", "#898781", "#fcfcfb"
C_FAV, C_POOR = "#2a78d6", "#eb6834"      # categorical slots 1 and 2
C_GOOD, C_BAD, C_WARN = "#0ca30c", "#d03b3b", "#eda100"
plt.rcParams.update({"figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
                     "axes.edgecolor": MUTED, "axes.labelcolor": INK2, "xtick.color": MUTED, "ytick.color": MUTED,
                     "text.color": INK, "font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.titlesize": 9.5, "axes.titleweight": "bold", "legend.frameon": False})


def km_panel(ax, groups, time, event, order, colors, title, xlabel="Days"):
    df = align_survival(groups, time, event)
    for g, c in zip(order, colors):
        sub = df[df.group == g]
        if sub.empty:
            continue
        kmf = KaplanMeierFitter().fit(sub.time, sub.event)
        sf = kmf.survival_function_.iloc[:, 0]
        ax.step(sf.index, sf.values, where="post", color=c, lw=2, label=f"{g} (n={len(sub)})")
        cens = sub[sub.event == 0].time
        ax.plot(cens, kmf.survival_function_at_times(cens).values, "|", color=c, ms=6, mew=1.2)
    p = logrank_p(df)
    ax.text(0.97, 0.95, f"log-rank p = {p:.3g}", transform=ax.transAxes, ha="right", va="top", color=INK2)
    ax.set_ylim(0, 1.03)
    ax.set_title(title, loc="left")
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Survival probability")
    ax.grid(axis="y", color="#e6e5e0", lw=0.6)
    ax.legend(loc="lower left")


def fig_km_dog():
    a = pd.read_csv(RESULTS / "cluster_assignments_dog.csv", index_col=0)
    a.index = a.index.astype(str)
    m = pd.read_pickle(W / "dog_meta.pkl")
    m.index = m.index.astype(str)
    fig, axes = plt.subplots(2, 2, figsize=(9, 6.6), constrained_layout=True)
    for i, s in enumerate(["GS-1", "GS-2"]):
        for j, (ep, col) in enumerate([("Disease-free interval", "dfi"), ("Overall survival", "os")]):
            km_panel(axes[i, j], a[f"{s}_rederived"], m[col], m[f"{col}_event"], ["FP", "PP"], [C_FAV, C_POOR],
                     f"DOG2, {s} clusters: {ep}")
            axes[i, j].set_xlim(0, 1600)
    fig.suptitle("Canine cohort (n = 186): published signatures, re-derived clusters", x=0.01, ha="left",
                 fontweight="bold")
    fig.savefig(FIGURES / "km_dog_gs.png", dpi=150)
    plt.close(fig)


def fig_km_target():
    a = pd.read_csv(RESULTS / "cluster_assignments_target.csv", index_col=0)
    m = pd.read_pickle(W / "hum_meta.pkl")
    fig, axes = plt.subplots(3, 4, figsize=(15, 10), constrained_layout=True)
    for i, (tag, lab) in enumerate([("hum", "all, n = 79"), ("hum_nonmet", "non-metastatic, n = 59"),
                                    ("hum_met", "metastatic, n = 20")]):
        for j, (s, ep, col) in enumerate([("GS-1", "PFS", "pfs"), ("GS-1", "OS", "os"),
                                          ("GS-2", "PFS", "pfs"), ("GS-2", "OS", "os")]):
            g = a[f"{s}_{tag}_rederived"].dropna()
            km_panel(axes[i, j], g, m.loc[g.index, col], m.loc[g.index, f"{col}_event"], ["FP", "PP"],
                     [C_FAV, C_POOR], f"TARGET {lab}: {s}, {ep}")
    fig.suptitle("Human TARGET-OS: canine GS-1 / GS-2 clusters, re-derived", x=0.01, ha="left", fontweight="bold")
    fig.savefig(FIGURES / "km_target_gs.png", dpi=130)
    plt.close(fig)


def fig_km_isa():
    a = pd.read_csv(RESULTS / "cluster_assignments_target_isa.csv", index_col=0)
    m = pd.read_pickle(W / "hum_meta.pkl")
    sigs = [c.rsplit("_", 1)[0] for c in a.columns if c.endswith("_hum")]
    fig, axes = plt.subplots(2, len(sigs), figsize=(3.6 * len(sigs), 6.6), constrained_layout=True, squeeze=False)
    for j, s in enumerate(sigs):
        for i, (tag, lab) in enumerate([("hum", "all"), ("hum_nonmet", "non-met")]):
            g = a[f"{s}_{tag}"].dropna()
            km_panel(axes[i, j], g, m.loc[g.index, "pfs"], m.loc[g.index, "pfs_event"], ["high", "low"],
                     [C_FAV, C_POOR], f"{s}, TARGET {lab}: PFS")
    fig.suptitle("Canine ISA signatures applied to TARGET-OS (k-means, k = 2)", x=0.01, ha="left",
                 fontweight="bold")
    fig.savefig(FIGURES / "km_target_isa.png", dpi=120)
    plt.close(fig)


def fig_heatmaps():
    gs = published_gs(D / "supp" / "42003_2023_5208_MOESM5_ESM.xlsx")
    omap = pd.read_pickle(W / "ortholog_map.pkl")
    da = pd.read_csv(RESULTS / "cluster_assignments_dog.csv", index_col=0)
    da.index = da.index.astype(str)
    ta = pd.read_csv(RESULTS / "cluster_assignments_target.csv", index_col=0)
    dt = pd.read_pickle(W / "dog_tpm_qn.pkl")
    dt.columns = dt.columns.astype(str)
    ht = pd.read_pickle(W / "hum_tpm_qn.pkl")
    panels = [("DOG2, GS-1", dt, gs["GS-1"], da["GS-1_rederived"]),
              ("TARGET, GS-1", ht, [h for h in map_signature(gs["GS-1-human"], omap).values() if h],
               ta["GS-1_hum_rederived"]),
              ("DOG2, ISA-3 genes", dt, PUBLISHED_ISA["ISA-3"], da["GS-1_rederived"]),
              ("TARGET, ISA-3 genes", ht, [h for h in map_signature(PUBLISHED_ISA["ISA-3"], omap).values() if h],
               ta["GS-1_hum_rederived"])]
    fig, axes = plt.subplots(1, 4, figsize=(16, 6), constrained_layout=True)
    for ax, (title, expr, genes, grp) in zip(axes, panels):
        z = log_zscore(expr)
        g, _ = signature_genes_present(genes, z.columns)
        zz = z[g].clip(-3, 3)
        order = grp.reindex(zz.index).sort_values().index
        gene_order = zz.loc[order].corrwith(pd.Series((grp[order] == "FP").astype(float).values, index=order))
        zz = zz.loc[order, gene_order.sort_values(ascending=False).index]
        im = ax.imshow(zz.T.values, aspect="auto", cmap="RdBu_r", vmin=-3, vmax=3, interpolation="nearest")
        n_fp = int((grp[order] == "FP").sum())
        ax.axvline(n_fp - 0.5, color=INK, lw=1)
        ax.set_yticks(range(len(zz.columns)))
        ax.set_yticklabels(zz.columns, fontsize=6)
        ax.set_xticks([n_fp / 2, n_fp + (len(order) - n_fp) / 2])
        ax.set_xticklabels([f"FP (n={n_fp})", f"PP (n={len(order) - n_fp})"])
        ax.set_title(title, loc="left")
    fig.colorbar(im, ax=axes, shrink=0.6, label="z-score of log2(TPM + 2)")
    fig.suptitle("Signature expression ordered by GS-1 cluster (samples) and correlation with FP (genes)",
                 x=0.01, ha="left", fontweight="bold")
    fig.savefig(FIGURES / "heatmaps_signatures.png", dpi=130)
    plt.close(fig)


def fig_deconv():
    dec = pd.read_csv(RESULTS / "deconvolution_fp_vs_pp.csv")
    da = pd.read_csv(RESULTS / "cluster_assignments_dog.csv", index_col=0)
    da.index = da.index.astype(str)
    ta = pd.read_csv(RESULTS / "cluster_assignments_target.csv", index_col=0)
    cells = [("TIL10", "T.cells.CD8", "CD8+ T cells (TIL10)"), ("TIL10", "Macrophages.M1", "M1 macrophages (TIL10)"),
             ("TIL10", "Macrophages.M2", "M2 macrophages (TIL10)"), ("MCP", "CD8 T cells", "CD8 T cells (MCP)"),
             ("MCP", "Monocytic lineage", "Monocytic lineage (MCP)")]
    fig, axes = plt.subplots(2, len(cells), figsize=(16, 6.4), constrained_layout=True)
    for i, (sp, grp, lab) in enumerate([("dog", da["GS-1_rederived"], "DOG2"),
                                        ("human", ta["GS-1_hum_rederived"], "TARGET")]):
        for j, (meth, cell, title) in enumerate(cells):
            est = pd.read_csv(RESULTS / f"deconv_{meth.lower()}_{sp}.csv.gz", index_col=0)
            est.columns = est.columns.astype(str)
            grp.index = grp.index.astype(str)
            v = [est.loc[cell, grp.index[grp == g]].values for g in ["FP", "PP"]]
            ax = axes[i, j]
            bp = ax.boxplot(v, widths=0.5, patch_artist=True, showfliers=False,
                            medianprops={"color": INK, "lw": 1.5})
            for patch, c in zip(bp["boxes"], [C_FAV, C_POOR]):
                patch.set_facecolor(c)
                patch.set_alpha(0.35)
                patch.set_edgecolor(c)
            rng = np.random.default_rng(0)
            for k, (vals, c) in enumerate(zip(v, [C_FAV, C_POOR]), start=1):
                ax.scatter(k + rng.uniform(-0.15, 0.15, len(vals)), vals, s=6, color=c, alpha=0.6, lw=0)
            r = dec[(dec.species == sp) & (dec.method == meth) & (dec.cell == cell)].iloc[0]
            ax.set_xticks([1, 2])
            ax.set_xticklabels(["FP", "PP"])
            ax.set_title(f"{lab}: {title}", loc="left", fontsize=8.5)
            ax.text(0.97, 0.97, f"BH q = {r.q_bh:.2g}", transform=ax.transAxes, ha="right", va="top", color=INK2)
    fig.suptitle("Open deconvolution by GS-1 cluster (MCP-counter scores; TIL10 constrained least-squares fractions)",
                 x=0.01, ha="left", fontweight="bold")
    fig.savefig(FIGURES / "deconvolution_fp_vs_pp.png", dpi=130)
    plt.close(fig)


def fig_summary():
    t = pd.read_csv(RESULTS / "replication_table.csv")
    t = t[t.paper_p.notna() & t.ours_p.notna()].reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(8.5, 0.28 * len(t) + 1.6), constrained_layout=True)
    y = np.arange(len(t))[::-1]
    lp = lambda p: -np.log10(np.clip(p, 1e-6, 1))  # noqa: E731
    ax.axvline(lp(0.05), color=MUTED, lw=1, ls="--")
    ax.text(lp(0.05), len(t) - 0.2, " p = 0.05", color=MUTED, va="bottom", fontsize=8)
    colors = {"yes": C_GOOD, "no": C_BAD, "partial": C_WARN}
    for yi, (_, r) in zip(y, t.iterrows()):
        ax.plot([lp(r.paper_p), lp(r.ours_p)], [yi, yi], color="#d6d5cf", lw=1.5, zorder=1)
        ax.scatter(lp(r.paper_p), yi, s=40, facecolor=SURF, edgecolor=INK2, lw=1.4, zorder=2, marker="o")
        ax.scatter(lp(r.ours_p), yi, s=46, color=colors.get(r.agree, MUTED), zorder=3,
                   marker={"yes": "o", "partial": "D", "no": "X"}.get(r.agree, "o"))
    ax.set_yticks(y)
    ax.set_yticklabels([f"C{int(r.claim)}  {r.test}" for _, r in t.iterrows()], fontsize=7.5)
    ax.set_xlabel("-log10 log-rank p   (open circle: paper; filled: this replication)")
    handles = [plt.Line2D([], [], marker=m, ls="", color=c, label=l, ms=7) for l, c, m in
               [("agrees (same call at p < 0.05)", C_GOOD, "o"), ("same direction, call differs", C_WARN, "D"),
                ("does not agree", C_BAD, "X")]]
    ax.legend(handles=handles, loc="lower right", fontsize=8)
    ax.set_title("Replication agreement, survival tests", loc="left")
    fig.savefig(FIGURES / "replication_summary.png", dpi=150)
    plt.close(fig)


def main() -> None:
    fig_km_dog()
    fig_km_target()
    fig_km_isa()
    fig_heatmaps()
    fig_deconv()
    fig_summary()
    print("figures:", sorted(p.name for p in FIGURES.glob("*.png")))


if __name__ == "__main__":
    main()
