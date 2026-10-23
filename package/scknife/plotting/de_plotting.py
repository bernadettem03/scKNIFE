import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import textwrap
from matplotlib.lines import Line2D
from scipy.stats import spearmanr
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist


def volcano_rank_score(group_df):
    """
    Composite ranking score for selecting top features after
    significance filtering.
    Score = |Cohen's d| × -log10(adjP)
    """
    d = group_df["cohens_d"].abs().max()
    p = max(group_df["adj_p_value"].min(), np.finfo(float).tiny)
    return d * (-np.log10(p))


def create_beeswarm(
    df,
    x="cohens_d",
    y="cell_type",
    cohen_threshold=0.2,
    pval_threshold=0.05,
    logfc_threshold=np.log2(1.1),
    show_n_significant=True,
    title=None,
    condition1_color="#8B0000",
    condition2_color="#000080",
    ns_color="#808080",
    strip_size=5,
    strip_alpha=0.7,
    strip_jitter=0.25,
    figsize=(12, 10),
    save=None,
    dpi=300,
    order=None,
):
    df = df.copy()

    group1 = df["group1"].iloc[0]
    group2 = df["group2"].iloc[0]

    df["significance"] = "NS"
    sig_mask = (
        (df["adj_p_value"] < pval_threshold)
        & (np.abs(df["cohens_d"]) >= cohen_threshold)
        & (np.abs(df["log2FC"]) >= logfc_threshold)
    )

    df.loc[sig_mask & (df["cohens_d"] < 0), "significance"] = group1
    df.loc[sig_mask & (df["cohens_d"] > 0), "significance"] = group2

    if order is None:
        order = (
            df[sig_mask]
            .groupby("cell_type")["cohens_d"]
            .apply(lambda x: np.mean(np.abs(x)))
            .sort_values(ascending=False)
            .index
        )
        if len(order) == 0:
            order = df["cell_type"].cat.categories if hasattr(df["cell_type"], "cat") else df["cell_type"].unique()

    df["cell_type"] = pd.Categorical(df["cell_type"], categories=order, ordered=True)

    counts = df[sig_mask]["cell_type"].value_counts()

    palette = {group1: condition1_color, group2: condition2_color, "NS": ns_color}

    fig, ax = plt.subplots(figsize=figsize)

    sns.stripplot(
        data=df,
        x=x,
        y=y,
        hue="significance",
        palette=palette,
        jitter=strip_jitter,
        size=strip_size,
        alpha=strip_alpha,
        ax=ax,
        legend=False,
    )

    ax.axvline(0, color="black", linestyle="--", lw=1)
    ax.axvline(cohen_threshold, color="black", linestyle=":", lw=1)
    ax.axvline(-cohen_threshold, color="black", linestyle=":", lw=1)

    if show_n_significant:
        xmax = ax.get_xlim()[1]
        for i, ct in enumerate(order):
            ax.text(
                xmax * 0.98,
                i,
                f"n={counts.get(ct, 0)}",
                ha="right",
                va="center",
                fontsize=11,
            )

    legend_elements = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor=condition1_color, markersize=10, label=f"FDR < {pval_threshold} ({group1})"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor=condition2_color, markersize=10, label=f"FDR < {pval_threshold} ({group2})"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor=ns_color, markersize=10, label=f"FDR ≥ {pval_threshold}"),
    ]
    ax.legend(handles=legend_elements, loc="lower right", frameon=False)

    ax.set_xlabel("Cohen's D", fontsize=14)
    ax.set_ylabel("")
    if title is None:
        title = f"{group1} vs {group2}"
    ax.set_title(title, fontsize=16)

    plt.tight_layout()

    if save:
        os.makedirs(os.path.dirname(save), exist_ok=True)
        base, _ = os.path.splitext(save)
        for ext in ['.pdf', '.png']:
            plt.savefig(base + ext, dpi=dpi, bbox_inches="tight")

    sig_df = df[df["significance"] != "NS"].set_index(["cell_type", "feature"])

    return fig, ax, sig_df


def create_volcano(
    de_results,
    effect_col="cohens_d",
    effect_threshold=0.2,
    padj_threshold=0.05,
    cell_type=None,
    group1=None,
    group2=None,
    effect_title="Cohen's d",
    wrapped_title_length=40,
    save=None,
    dpi=300,
    top_n_labels=5,
):
    df = de_results.copy()

    if cell_type is not None:
        df = df[df["cell_type"] == cell_type]
    if group1 is not None:
        df = df[df["group1"] == group1]
    if group2 is not None:
        df = df[df["group2"] == group2]

    df["neg_log_padj"] = -np.log10(df["adj_p_value"])
    df["neg_log_padj"] = df["neg_log_padj"].replace(np.inf, 300)

    sig = (df["adj_p_value"] < padj_threshold) & (df[effect_col].abs() >= effect_threshold)

    fig, ax = plt.subplots(figsize=(10, 8))

    ax.scatter(
        df.loc[~sig, effect_col],
        df.loc[~sig, "neg_log_padj"],
        color="lightgray",
        alpha=0.6,
        label="Not significant",
    )

    if sig.any():
        ax.scatter(
            df.loc[sig, effect_col],
            df.loc[sig, "neg_log_padj"],
            color="red",
            alpha=0.85,
            label="Significant",
        )

        df_sig = df[sig].copy()
        df_sig["rank_score"] = df_sig[effect_col].abs() * df_sig["neg_log_padj"]
        top_features = df_sig.nlargest(top_n_labels, "rank_score")

        for _, row in top_features.iterrows():
            ax.annotate(
                row["feature"],
                (row[effect_col], row["neg_log_padj"]),
                fontsize=8,
                alpha=0.9,
                ha="center",
                va="bottom",
            )

    ax.axhline(-np.log10(padj_threshold), color="black", linestyle="--", lw=1)
    ax.axvline(effect_threshold, color="black", linestyle="--", lw=1)
    ax.axvline(-effect_threshold, color="black", linestyle="--", lw=1)

    ax.set_xlabel(effect_title)
    ax.set_ylabel("-log10(FDR)")

    title_parts = []
    if cell_type:
        title_parts.append(cell_type)
    if group1 and group2:
        title_parts.append(f"{group1} vs {group2}")

    title = "Volcano Plot" if not title_parts else " – ".join(title_parts)
    ax.set_title("\n".join(textwrap.wrap(title, wrapped_title_length)))

    xmax = np.ceil(np.nanmax(np.abs(df[effect_col])))
    ax.set_xlim(-xmax, xmax)

    ax.legend()
    ax.grid(True, linestyle=":", alpha=0.5)

    plt.tight_layout()

    if save:
        os.makedirs(os.path.dirname(save), exist_ok=True)
        base, _ = os.path.splitext(save)
        for ext in ['.pdf', '.png']:
            plt.savefig(base + ext, dpi=dpi, bbox_inches="tight")

    if sig.any():
        return fig, ax, df.loc[sig].sort_values(effect_col).feature.unique().tolist()
    return fig, ax, []


def plot_top_features_heatmap(
    de_combined,
    top_n=50,
    top_n_per_group=None,
    figsize=(14, 10),
    save=None,
    dpi=300,
    adj_p_cutoff=0.05,
    cohen_cutoff=0.2,
    logfc_cutoff=np.log2(1.1),
    effect_col="log2FC",
    effect_label="log2FC",
    apply_logfc_filter=True,
    selection_mode="pvalue",
    symmetric_percentile=None,
    prefix_columns=True,
    pathway_wrap_width=58,
):
    df = de_combined.copy()

    df["label"] = df["expansion_group"] + "_" + df["cell_type"]
    required = {"feature", "cell_type", "adj_p_value", "cohens_d", effect_col}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Heatmap input lacks required columns: {sorted(missing)}")

    if apply_logfc_filter:
        significant = get_significant(
            df, adj_p_cutoff=adj_p_cutoff,
            cohen_cutoff=cohen_cutoff,
            logfc_cutoff=logfc_cutoff,
        )
    else:
        significant = df[
            (df["adj_p_value"] < adj_p_cutoff)
            & (df["cohens_d"].abs() >= cohen_cutoff)
        ].copy()

    if significant.empty:
        print("No significant features for heatmap.")
        return None, None

    if selection_mode == "balanced_signed_composite":
        if top_n_per_group is not None:
            raise ValueError("Balanced signed selection cannot be combined with top_n_per_group")
        ranked = significant.copy()
        ranked["composite_score"] = (
            ranked["cohens_d"].abs()
            * -np.log10(ranked["adj_p_value"].clip(lower=np.finfo(float).tiny))
        )
        representatives = (
            ranked.sort_values(
                ["composite_score", "feature", "cell_type"],
                ascending=[False, True, True],
                kind="mergesort",
            )
            .drop_duplicates("feature", keep="first")
        )
        n_negative = top_n // 2
        n_positive = top_n - n_negative
        positive = representatives[representatives["cohens_d"] > 0].head(n_positive)
        negative = representatives[representatives["cohens_d"] < 0].head(n_negative)
        selected = pd.concat([positive, negative], ignore_index=True)
        if len(selected) < top_n:
            remaining = representatives[~representatives["feature"].isin(selected["feature"])]
            selected = pd.concat(
                [selected, remaining.head(top_n - len(selected))], ignore_index=True
            )
        top_features = pd.Index(selected["feature"].head(top_n))
        if top_features.has_duplicates:
            raise RuntimeError("Balanced heatmap selection produced duplicate pathways")
    elif selection_mode != "pvalue":
        raise ValueError(f"Unknown heatmap selection_mode: {selection_mode}")
    elif top_n_per_group is not None:
        top_by_group = []
        for group in significant["expansion_group"].unique():
            group_sig = significant[significant["expansion_group"] == group]
            scores = (
                group_sig
                .groupby("feature", sort=False)
                .apply(volcano_rank_score)
            )
            top_features_group = (
                scores
                .nlargest(top_n_per_group)
                .index
            )
            top_by_group.append(top_features_group)
        top_features = pd.Index(
            pd.unique(
                np.concatenate(
                    [idx.to_numpy() for idx in top_by_group]
                )
            )
        )
    else:
        top_features = (
            significant.groupby("feature")["adj_p_value"]
            .min()
            .nsmallest(top_n)
            .index
        )

    column_key = "label" if prefix_columns else "cell_type"
    pivot = df[df["feature"].isin(top_features)].pivot_table(
        index="feature",
        columns=column_key,
        values=effect_col,
        aggfunc="first",
    ).fillna(0)

    if selection_mode == "balanced_signed_composite":
        if len(pivot) != len(top_features):
            raise RuntimeError("Not every selected pathway is represented in the heatmap matrix")
        if len(pivot) > 1:
            row_linkage = linkage(
                pdist(pivot.to_numpy(float)), method="average", optimal_ordering=True
            )
            dendro_info = dendrogram(row_linkage, no_plot=True)
            pivot = pivot.iloc[dendro_info["leaves"]]

        values = pivot.to_numpy(float)
        finite = np.abs(values[np.isfinite(values)])
        percentile = 98 if symmetric_percentile is None else symmetric_percentile
        vmax = float(np.percentile(finite, percentile))
        if not np.isfinite(vmax) or np.isclose(vmax, 0):
            vmax = float(np.nanmax(finite)) if len(finite) else 1.0
        vmax = vmax if vmax > 0 else 1.0

        fig = plt.figure(figsize=figsize)
        grid = fig.add_gridspec(
            1, 4, width_ratios=[1.1, 4.2, 5.0, 0.35], wspace=0.035
        )
        ax_dendrogram = fig.add_subplot(grid[0, 0])
        ax_heatmap = fig.add_subplot(grid[0, 1])
        ax_labels = fig.add_subplot(grid[0, 2], sharey=ax_heatmap)
        ax_colorbar = fig.add_subplot(grid[0, 3])

        if len(pivot) > 1:
            dendrogram(
                row_linkage,
                orientation="left",
                no_labels=True,
                color_threshold=0,
                above_threshold_color="#666666",
                link_color_func=lambda _: "#666666",
                ax=ax_dendrogram,
            )
        ax_dendrogram.invert_yaxis()
        ax_dendrogram.axis("off")

        sns.heatmap(
            pivot,
            cmap="RdBu_r",
            vmin=-vmax,
            vmax=vmax,
            center=0,
            xticklabels=True,
            yticklabels=False,
            cbar=True,
            cbar_ax=ax_colorbar,
            cbar_kws={"label": effect_label},
            ax=ax_heatmap,
        )
        ax_heatmap.set_xlabel("Cell type", fontsize=9)
        ax_heatmap.set_ylabel("")
        ax_heatmap.set_title(
            "Positive d: responder effect > non-responder effect\n"
            "Negative d: non-responder effect > responder effect",
            fontsize=8,
            pad=8,
        )
        ax_heatmap.tick_params(axis="x", labelsize=7, rotation=45)
        for tick in ax_heatmap.get_xticklabels():
            tick.set_ha("right")
            tick.set_rotation_mode("anchor")

        ax_labels.set_xlim(0, 1)
        ax_labels.set_ylim(ax_heatmap.get_ylim())
        ax_labels.set_xticks([])
        ax_labels.set_yticks([])
        for row, pathway in enumerate(pivot.index.astype(str)):
            ax_labels.text(
                0,
                row + 0.5,
                textwrap.fill(pathway, width=pathway_wrap_width),
                ha="left",
                va="center",
                fontsize=7,
                linespacing=0.9,
            )
        for spine in ax_labels.spines.values():
            spine.set_visible(False)
        ax_labels.set_title(
            "Metabolic task",
            loc="left",
            fontsize=9,
            fontweight="bold",
            pad=8,
        )
        ax_colorbar.tick_params(labelsize=7, length=2)
        ax_colorbar.yaxis.label.set_size(8)
        fig.subplots_adjust(left=0.025, right=0.985, bottom=0.12, top=0.975)

        if save:
            os.makedirs(os.path.dirname(save), exist_ok=True)
            base, _ = os.path.splitext(save)
            for ext in [".pdf", ".png"]:
                fig.savefig(base + ext, dpi=dpi, bbox_inches="tight", facecolor="white")
        return fig, ax_heatmap

    g = sns.clustermap(
        pivot,
        cmap="RdBu_r",
        center=0,
        figsize=figsize,
        row_cluster=True,
        col_cluster=False,
        xticklabels=True,
        yticklabels=True,
        dendrogram_ratio=(0.1, 0),
    )

    if save:
        os.makedirs(os.path.dirname(save), exist_ok=True)
        base, _ = os.path.splitext(save)
        for ext in ['.pdf', '.png']:
            g.savefig(base + ext, dpi=dpi, bbox_inches="tight")

    return g.fig, g.ax_heatmap


def plot_top_up_down_heatmap(
    de_combined,
    top_n=50,
    pval_cutoff=0.05,
    figsize=(16, 24),
    save=None,
    dpi=300,
):
    df = de_combined.copy()
    df["label"] = df["expansion_group"] + "_" + df["cell_type"]

    sig = df[df["adj_p_value"] < pval_cutoff]
    if sig.empty:
        print("No significant features for up/down heatmap.")
        return None, None

    feat_mean = sig.groupby("feature")["log2FC"].mean()

    top_on = feat_mean.nlargest(top_n).index
    top_pre = feat_mean.nsmallest(top_n).index
    top_features = list(top_on) + list(top_pre)

    pivot = df[df["feature"].isin(top_features)].pivot_table(
        index="feature", columns="label", values="log2FC", aggfunc="first"
    ).fillna(0).reindex(list(top_on) + list(top_pre))

    fig, ax = plt.subplots(figsize=figsize)
    sns.heatmap(
        pivot, cmap="RdBu_r", center=0,
        xticklabels=True, yticklabels=True,
        cbar_kws={"label": "log2FC"},
        ax=ax,
    )
    ax.set_ylabel("Metabolic tasks")
    ax.set_xlabel("Cell type")
    plt.tight_layout()
    if save:
        os.makedirs(os.path.dirname(save), exist_ok=True)
        base, _ = os.path.splitext(save)
        for ext in ['.pdf', '.png']:
            fig.savefig(base + ext, dpi=dpi, bbox_inches="tight")
    plt.close()
    return fig, ax


def plot_top_r_nr_heatmap(
    de_combined,
    top_n=50,
    pval_cutoff=0.05,
    figsize=(16, 24),
    save=None,
    dpi=300,
    row_cluster=True,
    col_cluster=True,
):
    df = de_combined.copy()
    df["label"] = df["expansion_group"] + "_" + df["cell_type"]

    sig = df[df["adj_p_value"] < pval_cutoff]
    if sig.empty:
        print("No significant features for R/NR heatmap.")
        return None, None

    feat_group = sig.groupby(["feature", "expansion_group"])["log2FC"].mean()

    top_r = feat_group.xs("Responders", level="expansion_group").nlargest(top_n).index
    top_nr = feat_group.xs("Non-Responders", level="expansion_group").nlargest(top_n).index
    top_features = list(top_r) + list(top_nr)

    pivot = df[df["feature"].isin(top_features)].pivot_table(
        index="feature", columns="label", values="log2FC", aggfunc="first"
    ).fillna(0).reindex(list(top_r) + list(top_nr))

    g = sns.clustermap(
        pivot, cmap="RdBu_r", center=0,
        xticklabels=True, yticklabels=True,
        cbar_kws={"label": "log2FC"},
        row_cluster=row_cluster,
        col_cluster=col_cluster,
        figsize=figsize,
    )
    g.ax_heatmap.set_ylabel("Metabolic tasks")
    g.ax_heatmap.set_xlabel("Cell type")
    if save:
        os.makedirs(os.path.dirname(save), exist_ok=True)
        base, _ = os.path.splitext(save)
        for ext in ['.pdf', '.png']:
            g.savefig(base + ext, dpi=dpi, bbox_inches="tight")
    if save:
        plt.close(g.fig)
    return g


def get_significant(df, adj_p_cutoff=0.05, cohen_cutoff=0.2, logfc_cutoff=np.log2(1.1)):
    mask = (
        (df["adj_p_value"] < adj_p_cutoff)
        & (df["cohens_d"].abs() >= cohen_cutoff)
        & (df["log2FC"].abs() >= logfc_cutoff)
    )
    return df[mask].copy()


def plot_logical_ranking(
    knife_df,
    sccellfie_df,
    top_n=10,
    title=None,
    save=None,
    dpi=300,
):
    knife = knife_df.copy()
    sccellfie = sccellfie_df.copy()

    knife["source"] = "KNIFE"
    sccellfie["source"] = "scCellFie"

    all_de = pd.concat([knife, sccellfie], ignore_index=True)

    logical_score = (
        knife.groupby(["cell_type", "feature"])["cohens_d"]
        .mean()
        .reset_index()
    )
    logical_score["score"] = logical_score["cohens_d"].abs()

    top_combos = (
        logical_score.groupby("feature")["score"]
        .max()
        .nlargest(top_n)
        .index
    )

    plot_df = all_de[all_de["feature"].isin(top_combos)].copy()
    plot_df["neg_log_padj"] = -np.log10(plot_df["adj_p_value"].clip(lower=1e-300))
    plot_df["size"] = plot_df["neg_log_padj"]

    fig, axes = plt.subplots(1, 2, figsize=(18, max(8, top_n * 0.6)), sharey=True)

    for ax_idx, (src_name, src_df) in enumerate([("KNIFE", knife), ("scCellFie", sccellfie)]):
        ax = axes[ax_idx]
        sub = plot_df[plot_df["source"] == src_name]

        scatter = ax.scatter(
            x=sub["cell_type"],
            y=sub["feature"],
            s=sub["size"] * 15,
            c=sub["log2FC"],
            cmap="RdBu_r",
            edgecolors="gray",
            linewidth=0.5,
            alpha=0.8,
        )

        ax.set_title(src_name, fontsize=14)
        ax.set_xlabel("Cell Type")
        ax.set_ylabel("Metabolic Task" if ax_idx == 0 else "")

        plt.colorbar(scatter, ax=ax, label="log2FC", shrink=0.6)

    fig.suptitle(title or "KNIFE vs scCellFie: Logical Ranking", fontsize=16)
    plt.tight_layout()

    if save:
        os.makedirs(os.path.dirname(save), exist_ok=True)
        base, _ = os.path.splitext(save)
        for ext in ['.pdf', '.png']:
            plt.savefig(base + ext, dpi=dpi, bbox_inches="tight")

    return fig, axes


def analyze_overlap(
    knife_sig,
    sccellfie_sig,
    knife_full=None,
    sccellfie_full=None,
    group_label="",
    save=None,
    dpi=300,
):
    cell_types = sorted(set(knife_sig["cell_type"]).union(sccellfie_sig["cell_type"]))

    rows = []
    scatter_data = []
    for ct in cell_types:
        k_set = set(knife_sig[knife_sig["cell_type"] == ct]["feature"])
        s_set = set(sccellfie_sig[sccellfie_sig["cell_type"] == ct]["feature"])

        overlap = k_set & s_set
        only_k = k_set - s_set
        only_s = s_set - k_set

        jaccard = len(overlap) / len(k_set | s_set) if len(k_set | s_set) > 0 else 0

        # Directional concordance among overlapping significant features
        k_sig_ct = knife_sig[knife_sig["cell_type"] == ct].set_index("feature")
        s_sig_ct = sccellfie_sig[sccellfie_sig["cell_type"] == ct].set_index("feature")
        shared = k_sig_ct.index.intersection(s_sig_ct.index)
        if len(shared) > 0:
            agree = (k_sig_ct.loc[shared, "cohens_d"] * s_sig_ct.loc[shared, "cohens_d"]) > 0
            concordance = agree.mean() * 100
        else:
            concordance = np.nan

        # Spearman correlation using full DE results
        rho = np.nan
        scf_full_ct = None
        k_full = knife_full[knife_full["cell_type"] == ct] if knife_full is not None else knife_sig[knife_sig["cell_type"] == ct]
        s_full = sccellfie_full[sccellfie_full["cell_type"] == ct] if sccellfie_full is not None else sccellfie_sig[sccellfie_sig["cell_type"] == ct]
        k_full = k_full[["feature", "cohens_d"]].drop_duplicates("feature")
        s_full = s_full[["feature", "cohens_d"]].drop_duplicates("feature")
        merged = k_full.merge(s_full, on="feature", suffixes=("_knife", "_sccellfie"))
        if len(merged) >= 5:
            rho, _ = spearmanr(merged["cohens_d_knife"], merged["cohens_d_sccellfie"], nan_policy="omit")
            if np.isnan(rho):
                rho = np.nan
            scatter_data.append({
                "cell_type": ct,
                "knife_d": merged["cohens_d_knife"].values,
                "sccellfie_d": merged["cohens_d_sccellfie"].values,
                "rho": rho,
            })

        rows.append({
            "cell_type": ct,
            "KNIFE_only": len(only_k),
            "scCellFie_only": len(only_s),
            "overlap": len(overlap),
            "jaccard": jaccard,
            "directional_concordance": concordance,
            "spearman_rho": rho,
        })

    summary = pd.DataFrame(rows)
    summary["total"] = summary["KNIFE_only"] + summary["overlap"] + summary["scCellFie_only"]
    summary["shared_pct"] = (summary["overlap"] / summary["total"] * 100).fillna(0)
    summary["knife_only_pct"] = (summary["KNIFE_only"] / summary["total"] * 100).fillna(0)
    summary["sccellfie_only_pct"] = (summary["scCellFie_only"] / summary["total"] * 100).fillna(0)

    fig, axes = plt.subplots(3, 2, figsize=(16, 16))
    ax1, ax2 = axes[0, 0], axes[0, 1]
    ax3, ax4 = axes[1, 0], axes[1, 1]
    ax5, ax6 = axes[2, 0], axes[2, 1]

    x = np.arange(len(summary))
    w = 0.25

    # ── Row 1 Col 1: Grouped bar ──
    ax1.bar(x - w, summary["KNIFE_only"], w, label="KNIFE only", color="#8B0000", alpha=0.7)
    ax1.bar(x, summary["overlap"], w, label="Overlap", color="#2E8B57", alpha=0.7)
    ax1.bar(x + w, summary["scCellFie_only"], w, label="scCellFie only", color="#000080", alpha=0.7)
    ax1.set_xticks(x)
    ax1.set_xticklabels(summary["cell_type"], rotation=45, ha="right")
    ax1.set_ylabel("Number of significant features")
    ax1.set_title(f"Overlap by Cell Type ({group_label})")
    ax1.legend()

    # ── Row 1 Col 2: Jaccard + Spearman ρ side-by-side ──
    x2 = np.arange(len(summary))
    w2 = 0.3
    ax2.bar(x2 - w2 / 2, summary["jaccard"], w2, label="Jaccard", color="#2E8B57", alpha=0.8)
    ax2.bar(x2 + w2 / 2, summary["spearman_rho"], w2, label="Spearman ρ", color="#FF8C00", alpha=0.8)
    ax2.set_xticks(x2)
    ax2.set_xticklabels(summary["cell_type"], rotation=45, ha="right")
    ax2.set_ylabel("Similarity")
    ax2.set_title(f"Similarity Metrics ({group_label})")
    ax2.set_ylim(-1, 1)
    ax2.axhline(0, color="gray", linestyle="-", lw=0.5)
    ax2.legend()

    # ── Row 2 Col 1: Stacked bar (absolute) ──
    ax3.bar(summary["cell_type"], summary["KNIFE_only"], label="KNIFE only", color="#8B0000", alpha=0.7)
    ax3.bar(summary["cell_type"], summary["overlap"], bottom=summary["KNIFE_only"], label="Overlap", color="#2E8B57", alpha=0.7)
    ax3.bar(summary["cell_type"], summary["scCellFie_only"],
            bottom=summary["KNIFE_only"] + summary["overlap"], label="scCellFie only", color="#000080", alpha=0.7)
    ax3.set_xticks(x)
    ax3.set_xticklabels(summary["cell_type"], rotation=45, ha="right")
    ax3.set_ylabel("Number of significant features")
    ax3.set_title(f"Stacked Overlap ({group_label})")
    ax3.legend()

    # ── Row 2 Col 2: Stacked bar (percentages) ──
    ax4.bar(summary["cell_type"], summary["knife_only_pct"], label="KNIFE only", color="#8B0000", alpha=0.7)
    ax4.bar(summary["cell_type"], summary["shared_pct"], bottom=summary["knife_only_pct"], label="Overlap", color="#2E8B57", alpha=0.7)
    ax4.bar(summary["cell_type"], summary["sccellfie_only_pct"],
            bottom=summary["knife_only_pct"] + summary["shared_pct"], label="scCellFie only", color="#000080", alpha=0.7)
    ax4.set_xticks(x)
    ax4.set_xticklabels(summary["cell_type"], rotation=45, ha="right")
    ax4.set_ylabel("Percentage (%)")
    ax4.set_title(f"Stacked Overlap (%) ({group_label})")
    ax4.legend()

    # ── Row 3: Scatter: KNIFE d vs scCellFie d per cell type ──
    colors_ct = plt.cm.tab10(np.linspace(0, 1, len(scatter_data)))
    for idx, d in enumerate(scatter_data):
        ax5.scatter(d["knife_d"], d["sccellfie_d"], s=15, alpha=0.5,
                    color=colors_ct[idx], label=f'{d["cell_type"]} (ρ={d["rho"]:.2f})')
    ax5.axhline(0, color="gray", linestyle="-", lw=0.5)
    ax5.axvline(0, color="gray", linestyle="-", lw=0.5)
    ax5.set_xlabel("KNIFE Cohen's d")
    ax5.set_ylabel("scCellFie Cohen's d")
    ax5.set_title(f"KNIFE vs scCellFie by Cell Type ({group_label})")
    ax5.legend(fontsize=8, loc="upper left")
    ax5.grid(True, linestyle=":", alpha=0.3)

    # Plot ρ bar on ax6 for clarity
    rho_vals = summary.dropna(subset=["spearman_rho"])
    if not rho_vals.empty:
        colors_rho = plt.cm.RdYlGn(rho_vals["spearman_rho"].clip(-1, 1) / 2 + 0.5)
        ax6.barh(rho_vals["cell_type"], rho_vals["spearman_rho"], color=colors_rho)
        ax6.set_xlim(-1, 1)
        ax6.axvline(0, color="gray", linestyle="-", lw=0.5)
        ax6.set_xlabel("Spearman ρ")
        ax6.set_title(f"Spearman Correlation ({group_label})")
    else:
        ax6.set_visible(False)

    plt.tight_layout()

    if save:
        os.makedirs(os.path.dirname(save), exist_ok=True)
        base, _ = os.path.splitext(save)
        for ext in ['.pdf', '.png']:
            plt.savefig(base + ext, dpi=dpi, bbox_inches="tight")

    return summary, fig, (ax1, ax2)


def plot_metabolic_program(de_df, group_label, top_n=10, save=None, title=None):
    de = de_df.copy()
    de["abs_d"] = de["cohens_d"].abs()
    de["neg_log_padj"] = -np.log10(de["adj_p_value"].clip(lower=1e-300))

    top_features = de.groupby("feature")["abs_d"].max().nlargest(top_n).index
    plot_df = de[de["feature"].isin(top_features)].copy()

    n_feat = len(top_features)
    fig, ax = plt.subplots(figsize=(14, max(6, n_feat * 0.45)))
    scatter = ax.scatter(
        x=plot_df["cell_type"], y=plot_df["feature"],
        s=plot_df["neg_log_padj"] * 15,
        c=plot_df["cohens_d"], cmap="RdBu_r",
        edgecolors="gray", linewidth=0.5, alpha=0.8,
    )
    ax.set_title(title or f"Top {top_n} tasks by |d_{group_label}|", fontsize=14)
    ax.set_xlabel("Cell Type")
    ax.set_ylabel("Metabolic Task")
    plt.colorbar(scatter, ax=ax, label="Cohen's d", shrink=0.6)
    plt.tight_layout()
    if save:
        os.makedirs(os.path.dirname(save), exist_ok=True)
        base, _ = os.path.splitext(save)
        for ext in ['.pdf', '.png']:
            plt.savefig(base + ext, dpi=300, bbox_inches="tight")
    if save:
        plt.close()

    print(f"\n  Top {top_n} features by |d_{group_label}|:")
    for i, ft in enumerate(top_features[:20]):
        ct_max = de[de["feature"] == ft].loc[de[de["feature"] == ft]["abs_d"].idxmax(), "cell_type"]
        d_max = de[de["feature"] == ft]["abs_d"].max()
        print(f"    {i+1:2d}. {ft:55s}  |d|={d_max:.3f}  (top cell: {ct_max})")
    return fig, ax


def plot_quadrant_metabolic_tasks(de_R, de_NR, cohen_threshold=0.3, pval_threshold=0.05, top_n=15, save=None, title=None):
    merged = de_R.merge(de_NR, on=["cell_type", "feature"], suffixes=("_R", "_NR"))

    def _max_abs_signed(series):
        idx = series.abs().idxmax()
        return series.loc[idx]

    aggr = merged.groupby("feature").agg(
        d_R=("cohens_d_R", _max_abs_signed),
        d_NR=("cohens_d_NR", _max_abs_signed),
        p_R=("adj_p_value_R", "min"),
        p_NR=("adj_p_value_NR", "min"),
        n_cell_types=("cell_type", "count"),
    ).reset_index()

    sig_R = (aggr["d_R"].abs() >= cohen_threshold) & (aggr["p_R"] < pval_threshold)
    sig_NR = (aggr["d_NR"].abs() >= cohen_threshold) & (aggr["p_NR"] < pval_threshold)
    conditions = [
        sig_R & sig_NR,
        sig_R & ~sig_NR,
        ~sig_R & sig_NR,
        ~sig_R & ~sig_NR,
    ]
    labels = ["Both sig", "R only", "NR only", "Neither"]
    aggr["sig_group"] = np.select(conditions, labels, default="Neither")
    aggr["divergence"] = (aggr["d_R"] - aggr["d_NR"]).abs()
    aggr["composite_score"] = aggr["divergence"] * (-np.log10(
        (aggr["p_R"] * aggr["p_NR"]).clip(lower=1e-300)
    ))
    top_features = aggr.nlargest(top_n, "composite_score")

    fig, ax = plt.subplots(figsize=(10, 10))
    colors = {"Both sig": "#d62728", "R only": "#ff7f0e", "NR only": "#1f77b4", "Neither": "#cccccc"}
    for group in labels:
        mask = aggr["sig_group"] == group
        if not mask.any():
            continue
        dot_sizes = -np.log10(
            (aggr.loc[mask, "p_R"] * aggr.loc[mask, "p_NR"]).clip(lower=1e-300)
        )
        dot_sizes = dot_sizes.clip(upper=10) * 15
        ax.scatter(
            aggr.loc[mask, "d_R"], aggr.loc[mask, "d_NR"],
            c=colors[group], label=f"{group} ({mask.sum()})",
            s=dot_sizes, alpha=0.6, edgecolors="none",
        )

    for _, row in top_features.iterrows():
        ax.annotate(
            row["feature"], (row["d_R"], row["d_NR"]),
            fontsize=6, alpha=0.8, ha="center", va="bottom",
            arrowprops=dict(arrowstyle="-", color="gray", lw=0.3),
        )

    lim = max(aggr["d_R"].abs().max(), aggr["d_NR"].abs().max()) * 1.15
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.axhline(0, color="black", linestyle="--", linewidth=0.8, alpha=0.5)
    ax.axvline(0, color="black", linestyle="--", linewidth=0.8, alpha=0.5)
    ax.plot([-lim, lim], [-lim, lim], color="gray", linestyle=":", linewidth=0.5, alpha=0.4)

    ax.set_xlabel("Cohen's d — Responders (Pre vs On)", fontsize=12)
    ax.set_ylabel("Cohen's d — Non-Responders (Pre vs On)", fontsize=12)
    if title is None:
        title = "Metabolic Task Quadrant Plot (scKNIFE)"
    ax.set_title(title, fontsize=14)
    ax.legend(fontsize=8, loc="upper left", title="Significance", title_fontsize=8)

    q1 = ((aggr["d_R"] >= 0) & (aggr["d_NR"] >= 0)).sum()
    q2 = ((aggr["d_R"] < 0) & (aggr["d_NR"] >= 0)).sum()
    q3 = ((aggr["d_R"] < 0) & (aggr["d_NR"] < 0)).sum()
    q4 = ((aggr["d_R"] >= 0) & (aggr["d_NR"] < 0)).sum()

    ax.text(lim*0.7, lim*0.88, f"Q1: BOTH ↑ (n={q1})", fontsize=9, ha="center",
            alpha=0.7, fontweight="bold")
    ax.text(lim*0.7, lim*0.80, "R increases, NR increases", fontsize=7, ha="center", alpha=0.5)
    ax.text(lim*0.7, lim*0.73, "Shared activation", fontsize=6, ha="center", alpha=0.4, style="italic")

    ax.text(-lim*0.3, lim*0.88, f"Q2: NR ↑ R ↓ (n={q2})", fontsize=9, ha="center",
            alpha=0.7, fontweight="bold", color="#2166AC")
    ax.text(-lim*0.3, lim*0.80, "NR-enriched", fontsize=7, ha="center", alpha=0.5)
    ax.text(-lim*0.3, lim*0.73, "Non-responder specific", fontsize=6, ha="center", alpha=0.4, style="italic")

    ax.text(-lim*0.7, -lim*0.73, f"Q3: BOTH ↓ (n={q3})", fontsize=9, ha="center",
            alpha=0.7, fontweight="bold")
    ax.text(-lim*0.7, -lim*0.80, "R decreases, NR decreases", fontsize=7, ha="center", alpha=0.5)
    ax.text(-lim*0.7, -lim*0.88, "Shared suppression", fontsize=6, ha="center", alpha=0.4, style="italic")

    ax.text(lim*0.3, -lim*0.73, f"Q4: R ↑ NR ↓ (n={q4})", fontsize=9, ha="center",
            alpha=0.7, fontweight="bold", color="#B2182B")
    ax.text(lim*0.3, -lim*0.80, "R-enriched ★", fontsize=7, ha="center", alpha=0.5)
    ax.text(lim*0.3, -lim*0.88, "Responder specific", fontsize=6, ha="center", alpha=0.4, style="italic")

    fig.text(0.5, -0.02,
             "Each point summarizes one metabolic task aggregated across cell types "
             "(max |Cohen's d|). Dot size encodes statistical confidence "
             "(larger = lower combined p-value).",
             ha='center', fontsize=9, style='italic')

    plt.tight_layout()
    if save:
        os.makedirs(os.path.dirname(save), exist_ok=True)
        base, _ = os.path.splitext(save)
        for ext in ['.pdf', '.png']:
            plt.savefig(base + ext, dpi=300, bbox_inches="tight")
    if save:
        plt.close()
        print(f"  Saved quadrant plot: {save}")
    return aggr


def plot_ranking_divergence(de_R, de_NR, scf_R, scf_NR, top_n=10, save=None):
    RENAME = {"cohens_d": "cohens_d_R", "log2FC": "log2FC_R", "adj_p_value": "adj_p_R"}
    de_R2 = de_R.rename(columns=RENAME)
    de_NR2 = de_NR.rename(columns={"cohens_d": "cohens_d_NR", "log2FC": "log2FC_NR", "adj_p_value": "adj_p_NR"})

    merged = de_R2.merge(de_NR2, on=["cell_type", "feature"])
    merged["diff_score"] = (merged["cohens_d_R"] - merged["cohens_d_NR"]).abs()

    top_features = merged.groupby("feature")["diff_score"].max().nlargest(top_n).index

    knife_both = pd.concat([de_R.assign(group="R"), de_NR.assign(group="NR")], ignore_index=True)
    knife_top = knife_both[knife_both["feature"].isin(top_features)].copy()
    knife_top["neg_log_padj"] = -np.log10(knife_top["adj_p_value"].clip(lower=1e-300))

    scf_both = pd.concat([scf_R.assign(group="R"), scf_NR.assign(group="NR")], ignore_index=True)
    scf_top = scf_both[scf_both["feature"].isin(top_features)].copy()
    scf_top["neg_log_padj"] = -np.log10(scf_top["adj_p_value"].clip(lower=1e-300))

    fig, axes = plt.subplots(1, 2, figsize=(18, max(8, top_n * 0.6)), sharey=True)
    for ax_idx, (src_name, src_df) in enumerate([("scKNIFE", knife_top), ("scCellFie", scf_top)]):
        ax = axes[ax_idx]
        scatter = ax.scatter(
            x=src_df["cell_type"], y=src_df["feature"],
            s=src_df["neg_log_padj"] * 15,
            c=src_df["cohens_d"], cmap="RdBu_r",
            edgecolors="gray", linewidth=0.5, alpha=0.8,
        )
        ax.set_title(f"{src_name} (color=d_R)", fontsize=14)
        ax.set_xlabel("Cell Type" if ax_idx == 1 else "")
        ax.set_ylabel("Metabolic Task" if ax_idx == 0 else "")
        plt.colorbar(scatter, ax=ax, label="Cohen's d", shrink=0.6)
    fig.suptitle(f"Top {top_n} tasks by |d_R - d_NR| (scKNIFE)", fontsize=16)
    plt.tight_layout()
    if save:
        os.makedirs(os.path.dirname(save), exist_ok=True)
        base, _ = os.path.splitext(save)
        for ext in ['.pdf', '.png']:
            plt.savefig(base + ext, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {save}")
    return merged


def plot_response_specific_ranking(
    knife_R, knife_NR, scf_R, scf_NR,
    top_n=20,
    adj_p_cutoff=0.05,
    cohen_cutoff=0.3,
    logfc_cutoff=0.25,
    figsize=(24, 12),
    save=None,
    dpi=300,
):
    """
    Plot response-specific metabolic task ranking comparing scKNIFE vs scCellFie.

    Composite ranking score identifies pathways that are simultaneously:
    - Biologically different between R and NR (divergence: |d_R - d_NR|)
    - Strongly altered by therapy (strength: |d_R| + |d_NR|)
    - Statistically well supported (significance: -log10(p_R × p_NR))

    Score = |d_R - d_NR| × (|d_R| + |d_NR|) × -log10(p_R × p_NR)

    Parameters
    ----------
    knife_R : pd.DataFrame
        scKNIFE DE results for Responders
    knife_NR : pd.DataFrame
        scKNIFE DE results for Non-Responders
    scf_R : pd.DataFrame
        scCellFie DE results for Responders
    scf_NR : pd.DataFrame
        scCellFie DE results for Non-Responders
    top_n : int
        Number of top features to display (default: 20)
    adj_p_cutoff : float
        Adjusted p-value threshold for significance filtering
    cohen_cutoff : float
        Minimum |Cohen's d| for significance filtering
    logfc_cutoff : float
        Minimum |log2FC| for significance filtering
    figsize : tuple
        Figure size
    save : str
        Path to save figure (without extension)
    dpi : int
        Resolution for saved figure

    Returns
    -------
    pd.DataFrame
        Merged data with response-specific scores
    """

    def _compute_score(merged):
        d_r = merged["cohens_d_R"]
        d_nr = merged["cohens_d_NR"]
        p_r = merged["adj_p_value_R"].clip(lower=1e-300)
        p_nr = merged["adj_p_value_NR"].clip(lower=1e-300)
        divergence = (d_r - d_nr).abs()
        strength = d_r.abs() + d_nr.abs()
        significance = -np.log10(p_r * p_nr)
        return divergence * strength * significance

    def _prepare_method(r_df, nr_df):
        r_sig = get_significant(r_df, adj_p_cutoff=adj_p_cutoff,
                                cohen_cutoff=cohen_cutoff, logfc_cutoff=logfc_cutoff)
        nr_sig = get_significant(nr_df, adj_p_cutoff=adj_p_cutoff,
                                 cohen_cutoff=cohen_cutoff, logfc_cutoff=logfc_cutoff)

        r_sel = r_sig[["cell_type", "feature", "cohens_d", "log2FC", "adj_p_value"]].copy()
        nr_sel = nr_sig[["cell_type", "feature", "cohens_d", "log2FC", "adj_p_value"]].copy()
        r_sel.columns = ["cell_type", "feature", "cohens_d_R", "log2FC_R", "adj_p_value_R"]
        nr_sel.columns = ["cell_type", "feature", "cohens_d_NR", "log2FC_NR", "adj_p_value_NR"]

        merged = r_sel.merge(nr_sel, on=["cell_type", "feature"])
        if merged.empty:
            return merged, r_sig, nr_sig

        merged["score"] = _compute_score(merged)
        return merged, r_sig, nr_sig

    methods = [
        {"name": "scKNIFE", "R": knife_R, "NR": knife_NR},
        {"name": "scCellFie", "R": scf_R, "NR": scf_NR},
    ]

    all_merged = []
    method_data = []

    for m in methods:
        merged, r_sig, nr_sig = _prepare_method(m["R"], m["NR"])
        if merged.empty:
            method_data.append({"name": m["name"], "merged": merged, "top_features": []})
            continue

        top_features = (
            merged.groupby("feature")["score"]
            .max()
            .nlargest(top_n)
            .index
        )

        r_plot = m["R"][m["R"]["feature"].isin(top_features)].copy()
        nr_plot = m["NR"][m["NR"]["feature"].isin(top_features)].copy()
        r_plot["group"] = "R"
        nr_plot["group"] = "NR"
        combined = pd.concat([r_plot, nr_plot], ignore_index=True)

        method_data.append({
            "name": m["name"],
            "merged": merged,
            "top_features": list(top_features),
            "combined": combined,
        })
        all_merged.append(merged)

    if all_merged:
        best_method_idx = max(range(len(all_merged)),
                              key=lambda i: all_merged[i]["score"].max())
        best_features = method_data[best_method_idx]["top_features"]
    else:
        best_features = []

    fig, axes = plt.subplots(1, 2, figsize=figsize)
    sc = None

    for idx, md in enumerate(method_data):
        ax = axes[idx]
        if not md["top_features"]:
            ax.set_visible(False)
            continue

        combined = md.get("combined")
        if combined is None:
            ax.set_visible(False)
            continue

        top_features = md["top_features"]
        celltypes = sorted(combined["cell_type"].unique())

        # Data density dot scaling
        n_features = len(top_features)
        n_celltypes = len(celltypes)
        cell_width = figsize[0] / max(n_celltypes, 1)
        cell_height = figsize[1] / max(n_features, 1)
        cell_area_points = (cell_width * 72) * (cell_height * 72)
        target_dot_size = cell_area_points * 0.15
        dot_scale = target_dot_size / 10

        for i, feat in enumerate(top_features):
            for j, ct in enumerate(celltypes):
                data = combined[(combined["feature"] == feat) & (combined["cell_type"] == ct)]
                for _, row in data.iterrows():
                    x_pos = j - 0.22 if row["group"] == "R" else j + 0.22
                    s_val = -np.log10(max(row["adj_p_value"], 1e-300))
                    if s_val > -np.log10(adj_p_cutoff):
                        s_val = min(s_val, 10)  # Cap at 10 (adj_p = 1e-10)
                        sc = ax.scatter(
                            x_pos, i,
                            s=s_val * dot_scale,
                            c=row["log2FC"],
                            cmap="RdBu_r", vmin=-2, vmax=2,
                            edgecolor="black", linewidth=0.6, alpha=0.85, zorder=3
                        )

        ax.set_xticks(range(len(celltypes)))
        ax.set_xticklabels(celltypes, rotation=45, ha="right", fontsize=11)
        ax.set_yticks(range(len(top_features)))
        ax.set_yticklabels([textwrap.fill(f, width=35) for f in top_features], fontsize=10)
        ax.set_title(md["name"], fontsize=16, fontweight="bold")
        ax.set_xlim(-0.5, len(celltypes) - 0.5)
        ax.set_ylim(-1, len(top_features))

        for j in range(len(celltypes)):
            ax.text(j - 0.22, -0.8, "R", ha="center", fontsize=9,
                    fontweight="bold", color="#2E8B57")
            ax.text(j + 0.22, -0.8, "NR", ha="center", fontsize=9,
                    fontweight="bold", color="#8B0000")
            ax.axvline(j + 0.5, color="grey", alpha=0.15)

    if sc is not None:
        cbar_ax = fig.add_axes([0.92, 0.3, 0.012, 0.4])
        fig.colorbar(sc, cax=cbar_ax, label="log2 Fold Change")

    fig.suptitle(
        "scKNIFE identifies stronger responder-specific metabolic programs "
        "than scCellFie following anti-PD1 therapy",
        fontsize=18, y=1.02
    )
    plt.tight_layout()

    if save:
        os.makedirs(os.path.dirname(save), exist_ok=True)
        base, _ = os.path.splitext(save)
        for ext in ['.pdf', '.png']:
            fig.savefig(base + ext, dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    if all_merged:
        return pd.concat(all_merged, ignore_index=True)
    return pd.DataFrame()
