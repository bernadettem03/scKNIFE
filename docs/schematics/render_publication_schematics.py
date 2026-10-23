#!/usr/bin/env python3
"""Render publication-quality scKNIFE schematics for documentation."""

from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Rectangle
import numpy as np

# ─── Figure parameters ──────────────────────────────────────────────────────
WORKFLOW_WIDTH_MM = 180
WORKFLOW_HEIGHT_MM = 145
CONNECTION_WIDTH_MM = 180
CONNECTION_HEIGHT_MM = 100
OUTPUT_DPI = 300
FONT_FAMILY = "Liberation Sans"
TITLE_PT = 10
STAGE_PT = 8.5
LABEL_PT = 7.5
SMALL_PT = 6.3
TINY_PT = 5.5
LINE_WIDTH = 1.15
ARROW_WIDTH = 1.25

INK = "#1F2937"
MUTED = "#596579"
OUTLINE = "#172033"
BLUE = "#2F80ED"
PURPLE = "#9B51E0"
GREEN = "#27AE60"
ORANGE = "#F2994A"
LIGHT_BLUE = "#EDF5FA"
LIGHT_PURPLE = "#F3EDF7"
LIGHT_GREEN = "#EEF7F1"
LIGHT_ORANGE = "#FFF5E9"
CELLTYPE_COLORS = ["#2F80ED", "#F2994A", "#27AE60", "#EB5757", "#9B51E0"]
SITE_COLORS = ["#377EB8", "#E41A1C", "#984EA3"]

ROOT = Path(__file__).resolve().parents[2]
OUTPUT_DIR = ROOT / "docs" / "assets"

plt.rcParams.update({
    "font.family": FONT_FAMILY,
    "font.size": LABEL_PT,
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
    "axes.grid": False,
})


def canvas(width_mm, height_mm):
    fig, ax = plt.subplots(figsize=(width_mm / 25.4, height_mm / 25.4))
    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    return fig, ax


def box(ax, x, y, w, h, edge=OUTLINE, fill="white", radius=0.018, lw=LINE_WIDTH):
    patch = FancyBboxPatch(
        (x, y), w, h,
        boxstyle=f"round,pad=0.008,rounding_size={radius}",
        linewidth=lw, edgecolor=edge, facecolor=fill,
    )
    ax.add_patch(patch)
    return patch


def arrow(ax, start, end, color=INK, connection="arc3,rad=0", style="-|>", ls="-"):
    patch = FancyArrowPatch(
        start, end, arrowstyle=style, mutation_scale=8,
        linewidth=ARROW_WIDTH, color=color, linestyle=ls,
        connectionstyle=connection, shrinkA=2, shrinkB=2,
    )
    ax.add_patch(patch)
    return patch


def text(ax, x, y, value, size=LABEL_PT, weight="normal", color=INK,
         ha="center", va="center", **kwargs):
    return ax.text(x, y, value, fontsize=size, fontweight=weight, color=color,
                   ha=ha, va=va, **kwargs)


def stage(ax, x0, x1, number, label):
    text(ax, x0, 0.965, f"{number}  {label.upper()}", size=STAGE_PT,
         weight="bold", color=INK, ha="left")
    ax.plot([x0, x1], [0.935, 0.935], color="#AEB7C4", lw=0.8)


def matrix(ax, x, y, w, h, color, rows=4, cols=5):
    gap = min(w / cols, h / rows) * 0.18
    cw = (w - gap * (cols - 1)) / cols
    ch = (h - gap * (rows - 1)) / rows
    values = np.array([[.85, .35, .65, .25, .75], [.3, .9, .45, .7, .4],
                       [.7, .25, .95, .5, .65], [.4, .8, .3, .7, .9]])
    for r in range(rows):
        for c in range(cols):
            ax.add_patch(Rectangle((x + c * (cw + gap), y + (rows - 1 - r) * (ch + gap)),
                                   cw, ch, facecolor=color, edgecolor="none",
                                   alpha=float(values[r % 4, c % 5])))


def network(ax, cx, cy, sx, sy, color=BLUE):
    points = np.array([[-.46, .05], [-.22, .42], [.12, .20], [.40, .44],
                       [.48, -.10], [.13, -.42], [-.28, -.34]])
    points[:, 0] = cx + points[:, 0] * sx
    points[:, 1] = cy + points[:, 1] * sy
    edges = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 0),
             (1, 6), (2, 5), (1, 5), (2, 4)]
    for i, j in edges:
        ax.plot(points[[i, j], 0], points[[i, j], 1], color=color, lw=0.75)
    for px, py in points:
        ax.add_patch(Circle((px, py), min(sx, sy) * .065, facecolor="white",
                            edgecolor=color, lw=0.85))


def save(fig, stem):
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    for suffix in ("png", "pdf", "svg"):
        kwargs = {"dpi": OUTPUT_DPI} if suffix == "png" else {}
        fig.savefig(OUTPUT_DIR / f"{stem}.{suffix}", facecolor="white", **kwargs)
    svg_path = OUTPUT_DIR / f"{stem}.svg"
    svg_path.write_text("\n".join(line.rstrip() for line in svg_path.read_text().splitlines()) + "\n")
    plt.close(fig)



def icon_scatter(ax, cx, cy, sx, sy, color=BLUE):
    pts = np.array([[-.42, -.12], [-.31, .20], [-.12, -.28], [.02, .10],
                    [.18, .30], [.30, -.08], [.42, .18], [.13, -.35]])
    for i, (px, py) in enumerate(pts):
        ax.add_patch(Circle((cx + px * sx, cy + py * sy), min(sx, sy) * .075,
                            facecolor=CELLTYPE_COLORS[i % len(CELLTYPE_COLORS)],
                            edgecolor="white", lw=.45))


def icon_sliders(ax, cx, cy, sx, sy, color=ORANGE):
    levels = [-.25, .15, -.05]
    for i, level in enumerate(levels):
        yy = cy + (i - 1) * sy * .28
        ax.plot([cx - sx * .42, cx + sx * .42], [yy, yy], color=OUTLINE, lw=.75)
        ax.add_patch(Circle((cx + level * sx, yy), min(sx, sy) * .09,
                            facecolor="white", edgecolor=color, lw=1.0))


def icon_cells(ax, cx, cy, sx, sy, color=BLUE):
    for i, (dx, dy) in enumerate([(-.28, .05), (0, .22), (.28, -.02), (0, -.25)]):
        c = CELLTYPE_COLORS[i]
        ax.add_patch(Circle((cx + dx * sx, cy + dy * sy), min(sx, sy) * .13,
                            facecolor=c + "22", edgecolor=c, lw=1.0))
        ax.add_patch(Circle((cx + dx * sx, cy + dy * sy), min(sx, sy) * .045,
                            facecolor=c, edgecolor="none"))


def icon_states(ax, cx, cy, sx, sy, color=BLUE):
    for i, dx in enumerate((-.28, 0, .28)):
        c = CELLTYPE_COLORS[i]
        ax.add_patch(Circle((cx + dx * sx, cy), min(sx, sy) * .15,
                            facecolor="white", edgecolor=c, lw=1.0))
        ax.plot([cx + dx * sx - sx * .07, cx + dx * sx + sx * .07],
                [cy, cy], color=c, lw=1.0 + .3 * i)


def icon_metabolic(ax, cx, cy, sx, sy, color=PURPLE):
    pts = [(cx - sx * .34, cy), (cx - sx * .12, cy + sy * .22),
           (cx + sx * .12, cy - sy * .15), (cx + sx * .35, cy + sy * .12)]
    for p1, p2 in zip(pts[:-1], pts[1:]):
        arrow(ax, p1, p2, color=color, style="-|>")
    for px, py in pts:
        ax.add_patch(Circle((px, py), min(sx, sy) * .075,
                            facecolor="white", edgecolor=color, lw=1.0))


def icon_signaling(ax, cx, cy, sx, sy, color=PURPLE):
    left, right = (cx - sx * .30, cy), (cx + sx * .30, cy)
    for px, c in [(left[0], BLUE), (right[0], PURPLE)]:
        ax.add_patch(Circle((px, cy), min(sx, sy) * .20,
                            facecolor="white", edgecolor=c, lw=1.1))
        ax.add_patch(Circle((px, cy), min(sx, sy) * .065,
                            facecolor=c, edgecolor="none"))
    arrow(ax, (left[0] + sx * .12, cy + sy * .06),
          (right[0] - sx * .12, cy + sy * .06), color=GREEN)
    arrow(ax, (right[0] - sx * .12, cy - sy * .08),
          (left[0] + sx * .12, cy - sy * .08), color=ORANGE)


def icon_enrichment(ax, cx, cy, sx, sy, color=GREEN):
    heights = [.30, .55, .82]
    for i, height in enumerate(heights):
        ax.add_patch(Rectangle((cx - sx * .34 + i * sx * .25, cy - sy * .32),
                               sx * .13, sy * height * .65,
                               facecolor=color, edgecolor="none", alpha=.55 + .18 * i))
    ax.plot([cx - sx * .40, cx + sx * .40], [cy - sy * .32, cy - sy * .32],
            color=OUTLINE, lw=.65)


def icon_metadata(ax, cx, cy, sx, sy, color=GREEN):
    for row, c in enumerate(SITE_COLORS):
        yy = cy + (.25 - row * .25) * sy
        ax.add_patch(Circle((cx - sx * .30, yy), min(sx, sy) * .065,
                            facecolor=c, edgecolor="none"))
        ax.plot([cx - sx * .15, cx + sx * .34], [yy, yy], color=OUTLINE, lw=.8)


def icon_card(ax, x, y, w, h, label, color, icon, fill="white", subtitle=None):
    box(ax, x, y, w, h, color, fill, radius=.011, lw=.9)
    icon(ax, x + w / 2, y + h * .62, w * .50, h * .34, color)
    text(ax, x + w / 2, y + h * .22, label, SMALL_PT, "bold")
    if subtitle:
        text(ax, x + w / 2, y + h * .075, subtitle, TINY_PT, color=MUTED)

def render_workflow():
    fig, ax = canvas(WORKFLOW_WIDTH_MM, WORKFLOW_HEIGHT_MM)

    def stage_panel(x, y, w, h, letter, title, side_fill=None):
        panel = FancyBboxPatch(
            (x, y), w, h,
            boxstyle="round,pad=0.012,rounding_size=.055",
            linewidth=1.8, edgecolor=INK, facecolor="white",
        )
        ax.add_patch(panel)
        if side_fill:
            side = FancyBboxPatch(
                (x + w * .84, y), w * .16, h,
                boxstyle="round,pad=0.012,rounding_size=.055",
                linewidth=0, facecolor=side_fill,
            )
            ax.add_patch(side)
        text(ax, x + .025, y + h - .035, f"{letter}  {title}", TITLE_PT,
             "bold", ha="left", va="top")

    # a. Inputs
    stage_panel(.025, .700, .950, .245, "a", "Inputs")
    matrix(ax, .105, .770, .105, .090, BLUE, rows=4, cols=5)
    text(ax, .158, .735, "Single-cell AnnData", LABEL_PT, "bold")
    text(ax, .158, .712, "cells × genes", TINY_PT, color=MUTED)

    text(ax, .640, .900, "Biological graph sources", LABEL_PT, "bold")
    source_labels = ["PPI", "Pathways", "Markers", "CellChat", "Metabolism", "Custom"]
    source_x = [.405, .500, .595, .690, .785, .880]
    for i, (cx, label) in enumerate(zip(source_x, source_labels)):
        if i in (1, 4):
            icon_metabolic(ax, cx, .810, .050, .060, ORANGE if i == 4 else GREEN)
        elif i == 2:
            icon_states(ax, cx, .810, .050, .060, PURPLE)
        elif i == 3:
            icon_signaling(ax, cx, .810, .050, .060, PURPLE)
        else:
            network(ax, cx, .810, .050, .060, BLUE)
        text(ax, cx, .745, label, TINY_PT, "bold")

    arrow(ax, (.500, .692), (.500, .650), INK)

    # b. Graph construction, tuning, and final factorization
    stage_panel(.025, .390, .950, .245, "b", "Model fitting", LIGHT_PURPLE)
    network(ax, .175, .500, .115, .115, BLUE)
    text(ax, .175, .435, "Unified graph G", LABEL_PT, "bold")
    text(ax, .175, .412, "selected or custom", TINY_PT, color=MUTED)

    icon_sliders(ax, .475, .505, .115, .100, ORANGE)
    text(ax, .475, .435, "Tune on subset", LABEL_PT, "bold")
    text(ax, .475, .412, "k · λ · βW · βH", TINY_PT, color=MUTED)

    matrix(ax, .735, .490, .050, .090, PURPLE, rows=4, cols=2)
    matrix(ax, .810, .515, .065, .045, PURPLE, rows=2, cols=5)
    text(ax, .760, .465, "W", LABEL_PT, "bold")
    text(ax, .842, .465, "H", LABEL_PT, "bold")
    text(ax, .800, .435, "scKNIFE", LABEL_PT, "bold")
    text(ax, .800, .412, "X ≈ Rg W H", TINY_PT, "bold", color=PURPLE)

    arrow(ax, (.245, .505), (.395, .505), BLUE)
    arrow(ax, (.555, .505), (.690, .505), PURPLE)
    text(ax, .325, .525, "fix G", TINY_PT, color=MUTED)
    text(ax, .622, .525, "lock parameters", TINY_PT, color=MUTED)

    arrow(ax, (.500, .382), (.500, .340), INK)

    # c. Six downstream capabilities with one icon scale and no nested boxes.
    stage_panel(.025, .055, .950, .270, "c", "Downstream analyses")
    tasks = [
        ("Clustering\n+ UMAP", BLUE, icon_scatter),
        ("Cell-type\nassignment", BLUE, icon_cells),
        ("Cell-state\nassignment", PURPLE, icon_states),
        ("Metabolic\ntasks", ORANGE, icon_metabolic),
        ("Cell-cell\nsignaling", PURPLE, icon_signaling),
        ("Pathway\nenrichment", GREEN, icon_enrichment),
    ]
    task_x = [.105, .260, .415, .570, .725, .880]
    for cx, (label, color, icon) in zip(task_x, tasks):
        icon(ax, cx, .190, .085, .085, color)
        text(ax, cx, .105, label, SMALL_PT, "bold")

    save(fig, "scknife-workflow")


def render_connection():
    fig, ax = canvas(CONNECTION_WIDTH_MM, CONNECTION_HEIGHT_MM)

    module_colors = [BLUE, PURPLE, GREEN, ORANGE, "#EB5757", "#56CCF2",
                     "#BB6BD9", "#6FCF97"]

    def program_modules(cx, cy, sx, sy):
        centers = np.array([[-.38, .24], [-.12, .30], [.15, .28], [.39, .20],
                            [-.34, -.20], [-.08, -.25], [.19, -.22], [.40, -.15]])
        offsets = np.array([[-.035, -.020], [.030, -.015], [.004, .035]])
        for i, center in enumerate(centers):
            points = np.column_stack((cx + (center[0] + offsets[:, 0]) * sx,
                                      cy + (center[1] + offsets[:, 1]) * sy))
            ax.plot(points[[0, 1, 2, 0], 0], points[[0, 1, 2, 0], 1],
                    color=module_colors[i], lw=.65, alpha=.65)
            for px, py in points:
                ax.add_patch(Circle((px, py), min(sx, sy) * .030,
                                    facecolor=module_colors[i], edgecolor="white",
                                    lw=.40))

    def mini_heatmap(x, y, w, h, palette):
        values = np.array([[.20, .75, .42, .92], [.84, .32, .68, .24],
                           [.48, .90, .28, .62], [.72, .46, .88, .35]])
        gap = .0028
        cw, ch = (w - 3 * gap) / 4, (h - 3 * gap) / 4
        for row in range(4):
            for col in range(4):
                ax.add_patch(Rectangle((x + col * (cw + gap),
                                        y + (3 - row) * (ch + gap)), cw, ch,
                                       facecolor=palette, edgecolor="none",
                                       alpha=.18 + .78 * values[row, col]))

    def count_matrix(x, y, w, h):
        sizes = np.array([[.25, .80, .45, .62], [.72, .30, .90, .40],
                          [.42, .68, .22, .82]])
        for row in range(3):
            for col in range(4):
                cx = x + (col + .5) * w / 4
                cy = y + (2.5 - row) * h / 3
                ax.add_patch(Circle((cx, cy), .003 + .007 * sizes[row, col],
                                    facecolor=CELLTYPE_COLORS[row],
                                    edgecolor="none", alpha=.82))

    # Quiet lane labels create hierarchy without enclosing the content.
    text(ax, .035, .935, "MODEL-DERIVED", SMALL_PT, "bold", color=PURPLE,
         ha="left")
    ax.plot([.035, .965], [.905, .905], color="#D8DEE8", lw=.8)
    text(ax, .035, .205, "METADATA GROUPING", SMALL_PT, "bold", color=GREEN,
         ha="left")
    ax.plot([.035, .965], [.180, .180], color="#D8DEE8", lw=.8)

    # W defines factor communities and their representative graph-node programs.
    matrix(ax, .040, .665, .090, .095, BLUE, rows=4, cols=5)
    text(ax, .085, .625, "W", TITLE_PT, "bold", color=BLUE)
    text(ax, .085, .597, "nodes × factors", TINY_PT, color=MUTED)

    program_modules(.245, .712, .200, .170)
    text(ax, .245, .625, "Factor communities", LABEL_PT, "bold")
    text(ax, .245, .597, "CC1–CC8 define node programs", TINY_PT, color=MUTED)

    matrix(ax, .390, .665, .090, .095, PURPLE, rows=4, cols=5)
    text(ax, .435, .625, "P", TITLE_PT, "bold", color=PURPLE)
    text(ax, .435, .597, "programs × factors", TINY_PT, color=MUTED)

    arrow(ax, (.137, .712), (.155, .712), BLUE)
    arrow(ax, (.335, .712), (.377, .712), PURPLE)

    # H carries the same factor axis for cells and meets P only at reconstruction.
    matrix(ax, .390, .365, .090, .095, PURPLE, rows=4, cols=5)
    text(ax, .435, .325, "H", TITLE_PT, "bold", color=PURPLE)
    text(ax, .435, .297, "cells × factors", TINY_PT, color=MUTED)

    box(ax, .555, .455, .155, .265, PURPLE, "#FBF8FD", radius=.015, lw=1.25)
    text(ax, .632, .676, "A = H × Pᵀ", TITLE_PT, "bold", color=PURPLE)
    matrix(ax, .582, .545, .101, .083, PURPLE, rows=4, cols=5)
    text(ax, .632, .505, "Per-cell program activity", SMALL_PT, "bold")
    text(ax, .632, .477, "cells × programs", TINY_PT, color=MUTED)

    arrow(ax, (.487, .712), (.547, .635), PURPLE)
    arrow(ax, (.487, .412), (.547, .525), PURPLE)

    # Analytical outputs are miniature data views, not generic destination cards.
    output_specs = [
        (.705, "NG × cell type", "mean activity A", PURPLE, "heat"),
        (.465, "CC × cell type", "mean H by community", BLUE, "heat"),
        (.265, "Cell type × patient", "observed cell counts", GREEN, "count"),
    ]
    for y, title_value, subtitle, color, kind in output_specs:
        ax.add_patch(Rectangle((.790, y - .060), .175, .145, facecolor="white",
                               edgecolor="#C9D1DC", lw=.8))
        ax.add_patch(Rectangle((.790, y - .060), .007, .145,
                               facecolor=color, edgecolor="none"))
        if kind == "heat":
            mini_heatmap(.810, y - .032, .052, .085, color)
        else:
            count_matrix(.807, y - .035, .058, .088)
        text(ax, .878, y + .025, title_value, SMALL_PT, "bold", ha="left")
        text(ax, .878, y - .013, subtitle, TINY_PT, color=MUTED, ha="left")

    # Solid routes carry model-derived values to the first two summaries.
    arrow(ax, (.718, .590), (.782, .705), PURPLE)
    ax.plot([.487, .745, .745], [.382, .382, .465], color=BLUE,
            lw=ARROW_WIDTH)
    arrow(ax, (.745, .465), (.782, .465), BLUE)

    # Metadata stays in a separate lower lane and joins outputs only for grouping.
    icon_cells(ax, .095, .105, .080, .075, GREEN)
    text(ax, .150, .125, "Cell identity", SMALL_PT, "bold", ha="left")
    text(ax, .150, .092, "join by cell ID", TINY_PT, color=MUTED, ha="left")
    arrow(ax, (.225, .108), (.270, .108), GREEN, ls="--")

    metadata_items = [("cell type", BLUE), ("H cluster", PURPLE),
                      ("patient", GREEN), ("site", ORANGE)]
    for i, (label, color) in enumerate(metadata_items):
        x = .285 + i * .095
        ax.add_patch(Circle((x, .108), .007, facecolor=color, edgecolor="none"))
        text(ax, x + .012, .108, label, TINY_PT, "bold", ha="left")

    ax.plot([.665, .985], [.108, .108], color=GREEN, lw=1.0, ls="--")
    ax.plot([.985, .985], [.108, .760], color=GREEN, lw=1.0, ls="--")
    for y in (.690, .450, .250):
        arrow(ax, (.985, y), (.968, y), GREEN, ls="--")

    # Compact legend reinforces the model/metadata boundary.
    ax.plot([.615, .650], [.035, .035], color=PURPLE, lw=1.4)
    text(ax, .657, .035, "model transformation", TINY_PT, color=MUTED, ha="left")
    ax.plot([.795, .830], [.035, .035], color=GREEN, lw=1.2, ls="--")
    text(ax, .837, .035, "metadata grouping", TINY_PT, color=MUTED, ha="left")

    save(fig, "scknife-ng-cc-metadata")


if __name__ == "__main__":
    render_workflow()
    render_connection()
    print(f"Wrote publication schematics to {OUTPUT_DIR}")
