#!/usr/bin/env python3
"""
Replicate the legacy Step 11 activity-based ligand-receptor communication model.

This is not the CellChat R probability model. It combines Scknife activities
with CellChatDB ligand/receptor definitions and the exact legacy weighting
semantics, including fixed observed specificity during label permutations.
"""

from __future__ import annotations

import argparse
import gc
import hashlib
import json
import math
import time
from pathlib import Path

import anndata as ad
import networkx as nx
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from statsmodels.stats.multitest import multipletests

from .common import (
    ACTIVITY_PATH,
    CELLCHAT_DB_PATH,
    DATA_DIR,
    DB_ANNOTATIONS,
    FDR_THRESHOLD,
    GROUP_COLUMN,
    H5AD_PATH,
    METASTATIC_COLUMN,
    N_PERMUTATIONS,
    PATIENT_COLUMN,
    RANDOM_SEED,
    SITE_COLUMN,
    atomic_npz,
    ensure_directories,
    software_versions,
    write_json,
)


COFACTOR_MULTIPLIERS = {
    "agonist": 1.5,
    "antagonist": 0.5,
    "co_A_receptor": 1.3,
    "co_I_receptor": 0.7,
}


def parse_gene_list(value) -> list[str]:
    if pd.isna(value):
        return []
    return [item.strip() for item in str(value).split(",") if item.strip()]


def load_interactions() -> tuple[pd.DataFrame, list[str]]:
    db = pd.read_csv(CELLCHAT_DB_PATH)
    required = {
        "ligand.symbol", "receptor.symbol", "annotation", "evidence",
        "agonist", "antagonist", "co_A_receptor", "co_I_receptor",
    }
    missing = sorted(required - set(db.columns))
    if missing:
        raise ValueError(f"CellChatDB is missing required columns: {missing}")

    db = db.loc[db["annotation"].isin(DB_ANNOTATIONS)].copy()
    db["ligand"] = db["ligand.symbol"].astype(str).str.strip()
    db["receptors"] = db["receptor.symbol"].map(parse_gene_list)
    for column in COFACTOR_MULTIPLIERS:
        db[f"{column}_genes"] = db[column].map(parse_gene_list)
    db["confidence"] = db["evidence"].map(
        lambda value: 1.0
        if pd.notna(value) and "KEGG" in str(value)
        else 0.8 if pd.notna(value) else 0.5
    )

    genes = set(db["ligand"])
    for values in db["receptors"]:
        genes.update(values)
    for column in COFACTOR_MULTIPLIERS:
        for values in db[f"{column}_genes"]:
            genes.update(values)
    genes.discard("")
    return db, sorted(genes)


def activity_schema() -> tuple[pq.ParquetFile, list[str], list[str]]:
    parquet = pq.ParquetFile(ACTIVITY_PATH)
    columns = list(parquet.schema.names)
    index_column = "__index_level_0__"
    if index_column not in columns:
        raise ValueError(f"{ACTIVITY_PATH} has no pandas index column")
    cell_columns = [column for column in columns if column != index_column]
    node_names = [
        str(value)
        for value in pq.read_table(ACTIVITY_PATH, columns=[index_column])
        .column(0)
        .to_pylist()
    ]
    if len(node_names) != parquet.metadata.num_rows:
        raise ValueError("Activity node index length does not match parquet rows")
    if len(node_names) != len(set(node_names)):
        raise ValueError("Activity node names are not unique")
    return parquet, node_names, cell_columns


def load_metadata(cell_columns: list[str]) -> tuple[pd.DataFrame, list[str]]:
    adata = ad.read_h5ad(H5AD_PATH, backed="r")
    required = {GROUP_COLUMN, METASTATIC_COLUMN, SITE_COLUMN, PATIENT_COLUMN}
    missing = sorted(required - set(adata.obs.columns))
    if missing:
        adata.file.close()
        raise ValueError(f"H5AD metadata is missing columns: {missing}")

    metadata = adata.obs.loc[:, sorted(required)].copy()
    adata.file.close()
    if metadata.index.has_duplicates:
        raise ValueError("H5AD cell barcodes are not unique")
    if len(cell_columns) != len(set(cell_columns)):
        raise ValueError("Activity parquet cell columns are not unique")

    missing_in_h5ad = sorted(set(cell_columns) - set(metadata.index))
    missing_in_activity = sorted(set(metadata.index) - set(cell_columns))
    if missing_in_h5ad or missing_in_activity:
        raise ValueError(
            "Barcode mismatch: "
            f"{len(missing_in_h5ad)} activity-only and "
            f"{len(missing_in_activity)} H5AD-only cells"
        )

    metadata = metadata.loc[cell_columns].copy()
    if metadata[GROUP_COLUMN].isna().any():
        raise ValueError(f"{GROUP_COLUMN} contains missing values")
    metadata[GROUP_COLUMN] = metadata[GROUP_COLUMN].astype(str)
    group_order = sorted(metadata[GROUP_COLUMN].unique().tolist())
    if len(group_order) < 2:
        raise ValueError(f"{GROUP_COLUMN} must contain at least two assigned groups")
    return metadata, group_order


def stream_activity_subset(
    wanted_genes: list[str],
    node_names: list[str],
    cell_columns: list[str],
    chunk_cells: int,
) -> tuple[np.ndarray, list[str], dict]:
    node_to_row = {node: row for row, node in enumerate(node_names)}
    available = [gene for gene in wanted_genes if gene in node_to_row]
    row_positions = np.array([node_to_row[gene] for gene in available], dtype=np.int64)
    activity = np.empty((len(cell_columns), len(available)), dtype=np.float32)
    read_start = time.time()

    for start in range(0, len(cell_columns), chunk_cells):
        stop = min(start + chunk_cells, len(cell_columns))
        columns = cell_columns[start:stop]
        table = pq.read_table(ACTIVITY_PATH, columns=columns)
        block = table.to_pandas().to_numpy(dtype=np.float32, copy=False)
        activity[start:stop, :] = block[row_positions, :].T
        del table, block
        gc.collect()
        print(f"  activity columns {stop:,}/{len(cell_columns):,}", flush=True)

    if not np.isfinite(activity).all():
        raise ValueError("Selected activity matrix contains non-finite values")
    if np.min(activity) < 0:
        raise ValueError("Selected activity matrix contains negative values")
    stats = {
        "wanted_nodes": len(wanted_genes),
        "available_nodes": len(available),
        "missing_nodes": sorted(set(wanted_genes) - set(available)),
        "minimum": float(activity.min()),
        "maximum": float(activity.max()),
        "elapsed_seconds": time.time() - read_start,
    }
    return activity, available, stats


def group_means(
    activity: np.ndarray,
    labels: np.ndarray,
    n_groups: int,
) -> tuple[np.ndarray, np.ndarray]:
    means = np.zeros((activity.shape[1], n_groups), dtype=np.float64)
    counts = np.bincount(labels, minlength=n_groups).astype(np.int64)
    if np.any(counts == 0):
        raise ValueError(f"At least one H cluster has no cells: {counts.tolist()}")
    for group in range(n_groups):
        means[:, group] = activity[labels == group].mean(axis=0, dtype=np.float64)
    return means, counts


def build_pairs(
    db: pd.DataFrame,
    available_genes: list[str],
    observed_means: np.ndarray,
) -> pd.DataFrame:
    gene_to_index = {gene: idx for idx, gene in enumerate(available_genes)}
    records = []
    global_gene_activity = observed_means.mean(axis=1)

    for db_index, row in db.iterrows():
        ligand = row["ligand"]
        if ligand not in gene_to_index:
            continue

        cofactor_weight = 1.0
        available_cofactors = []
        for column, multiplier in COFACTOR_MULTIPLIERS.items():
            genes = [
                gene for gene in row[f"{column}_genes"]
                if gene in gene_to_index
            ]
            if genes:
                activity_value = max(
                    global_gene_activity[gene_to_index[gene]] for gene in genes
                )
                if activity_value > 0:
                    cofactor_weight *= multiplier
                available_cofactors.extend(genes)

        for receptor in row["receptors"]:
            if receptor not in gene_to_index:
                continue
            records.append({
                "db_row": int(db_index),
                "interaction_name": row.get("interaction_name", ""),
                "pathway_name": row.get("pathway_name", ""),
                "annotation": row["annotation"],
                "evidence": row["evidence"],
                "ligand": ligand,
                "receptor": receptor,
                "ligand_index": gene_to_index[ligand],
                "receptor_index": gene_to_index[receptor],
                "confidence": float(row["confidence"]),
                "cofactor_weight": float(cofactor_weight),
                "available_cofactors": ",".join(sorted(set(available_cofactors))),
            })
    pairs = pd.DataFrame(records)
    if pairs.empty:
        raise ValueError("No valid ligand-receptor components were found")
    return pairs


def communication_matrix(
    means: np.ndarray,
    specificity: np.ndarray,
    pairs: pd.DataFrame,
) -> np.ndarray:
    ligand_idx = pairs["ligand_index"].to_numpy(dtype=np.int64)
    receptor_idx = pairs["receptor_index"].to_numpy(dtype=np.int64)
    pair_weight = (
        pairs["confidence"].to_numpy(dtype=np.float64)
        * pairs["cofactor_weight"].to_numpy(dtype=np.float64)
    )
    left = means[ligand_idx, :] * specificity[ligand_idx, :]
    right = means[receptor_idx, :] * specificity[receptor_idx, :]
    return left.T @ (right * pair_weight[:, None])


def literal_legacy_matrix(
    means: np.ndarray,
    specificity: np.ndarray,
    pairs: pd.DataFrame,
) -> np.ndarray:
    result = np.zeros((means.shape[1], means.shape[1]), dtype=np.float64)
    for row in pairs.itertuples(index=False):
        for source in range(means.shape[1]):
            for target in range(means.shape[1]):
                result[source, target] += (
                    means[row.ligand_index, source]
                    * means[row.receptor_index, target]
                    * specificity[row.ligand_index, source]
                    * specificity[row.receptor_index, target]
                    * row.confidence
                    * row.cofactor_weight
                )
    return result


def validate_vectorization(
    means: np.ndarray,
    specificity: np.ndarray,
    pairs: pd.DataFrame,
) -> float:
    test_pairs = pairs.head(min(20, len(pairs))).copy()
    fast = communication_matrix(means, specificity, test_pairs)
    literal = literal_legacy_matrix(means, specificity, test_pairs)
    error = float(np.max(np.abs(fast - literal)))
    if not np.allclose(fast, literal, rtol=1e-12, atol=1e-12):
        raise AssertionError(f"Vectorized score differs from legacy loops: {error}")
    return error


def interaction_details(
    means: np.ndarray,
    specificity: np.ndarray,
    pairs: pd.DataFrame,
    groups: list[str],
) -> pd.DataFrame:
    n_groups = len(groups)
    ligand_idx = pairs["ligand_index"].to_numpy(dtype=np.int64)
    receptor_idx = pairs["receptor_index"].to_numpy(dtype=np.int64)
    left = means[ligand_idx] * specificity[ligand_idx]
    right = means[receptor_idx] * specificity[receptor_idx]
    scores = (
        left[:, :, None]
        * right[:, None, :]
        * pairs["confidence"].to_numpy(dtype=np.float64)[:, None, None]
        * pairs["cofactor_weight"].to_numpy(dtype=np.float64)[:, None, None]
    )

    pair_index = np.repeat(np.arange(len(pairs)), n_groups * n_groups)
    source_index = np.tile(
        np.repeat(np.arange(n_groups), n_groups), len(pairs)
    )
    target_index = np.tile(np.arange(n_groups), len(pairs) * n_groups)
    details = pairs.iloc[pair_index].reset_index(drop=True).drop(
        columns=["ligand_index", "receptor_index"]
    )
    details.insert(0, "source", np.asarray(groups)[source_index])
    details.insert(1, "target", np.asarray(groups)[target_index])
    details["ligand_activity"] = means[ligand_idx[pair_index], source_index]
    details["receptor_activity"] = means[
        receptor_idx[pair_index], target_index
    ]
    details["specificity_weight"] = (
        specificity[ligand_idx[pair_index], source_index]
        * specificity[receptor_idx[pair_index], target_index]
    )
    details["weighted_score"] = scores.reshape(-1)
    return details.loc[details["weighted_score"] > 0].reset_index(drop=True)


def network_centrality(comm: pd.DataFrame) -> pd.DataFrame:
    graph = nx.DiGraph()
    for source in comm.index:
        for target in comm.columns:
            value = float(comm.loc[source, target])
            if value > 0:
                graph.add_edge(source, target, weight=value)
    out_degree = nx.out_degree_centrality(graph)
    in_degree = nx.in_degree_centrality(graph)
    betweenness = nx.betweenness_centrality(graph, weight="weight")
    try:
        information = nx.current_flow_betweenness_centrality(
            graph.to_undirected(), weight="weight"
        )
    except Exception:
        information = {node: 0.0 for node in graph.nodes}
    return pd.DataFrame({
        GROUP_COLUMN: list(graph.nodes),
        "out_degree": [out_degree.get(node, 0.0) for node in graph.nodes],
        "in_degree": [in_degree.get(node, 0.0) for node in graph.nodes],
        "betweenness": [betweenness.get(node, 0.0) for node in graph.nodes],
        "information": [information.get(node, 0.0) for node in graph.nodes],
    })


def run_permutations(
    activity: np.ndarray,
    labels: np.ndarray,
    observed: np.ndarray,
    specificity: np.ndarray,
    pairs: pd.DataFrame,
    n_groups: int,
    n_permutations: int,
    checkpoint_every: int,
    force: bool,
    label_signature: str,
) -> tuple[np.ndarray, int]:
    if n_permutations < 1 or checkpoint_every < 1:
        raise ValueError("Permutation counts and checkpoint interval must be positive")
    fingerprint = hashlib.sha256(b"scknife-permutation-v2")
    fingerprint.update(str((RANDOM_SEED, n_groups, label_signature)).encode())
    for values in (activity, labels, observed, specificity):
        values = np.ascontiguousarray(values)
        fingerprint.update(str((values.shape, values.dtype.str)).encode())
        fingerprint.update(values.tobytes())
    fingerprint.update(pairs.to_json(orient="split", double_precision=15).encode())
    input_signature = fingerprint.hexdigest()
    checkpoint = DATA_DIR / "permutation_checkpoint.npz"
    rng = np.random.default_rng(RANDOM_SEED)
    exceed = np.zeros_like(observed, dtype=np.int64)
    completed = 0

    if checkpoint.exists() and not force:
        saved = np.load(checkpoint, allow_pickle=False)
        saved_signature = str(saved["input_signature"]) if "input_signature" in saved.files else None
        compatible = (
            saved_signature == input_signature
            and saved["exceed"].shape == observed.shape
            and int(saved["completed"]) <= n_permutations
        )
        if compatible:
            completed = int(saved["completed"])
            exceed = saved["exceed"].astype(np.int64)
            rng.bit_generator.state = json.loads(str(saved["rng_state"]))
            print(f"Resuming compatible permutations at {completed}/{n_permutations}")
        else:
            print("Existing permutation checkpoint is incompatible with current inputs; restarting at 0")
    elif checkpoint.exists() and force:
        checkpoint.unlink()

    start = time.time()
    for permutation in range(completed, n_permutations):
        shuffled = rng.permutation(labels)
        null_means, _ = group_means(activity, shuffled, n_groups)
        null_comm = communication_matrix(null_means, specificity, pairs)
        exceed += null_comm >= observed
        completed = permutation + 1

        if completed % checkpoint_every == 0 or completed == n_permutations:
            atomic_npz(
                checkpoint,
                completed=np.array(completed, dtype=np.int64),
                exceed=exceed,
                rng_state=np.array(json.dumps(rng.bit_generator.state)),
                label_signature=np.array(label_signature),
                input_signature=np.array(input_signature),
            )
            elapsed = time.time() - start
            rate = completed / max(elapsed, 1e-9)
            print(
                f"  permutations {completed}/{n_permutations} "
                f"({rate:.2f}/s since process start)",
                flush=True,
            )
    return (exceed + 1) / float(n_permutations + 1), completed


def save_pathway_outputs(
    details: pd.DataFrame,
    groups: list[str],
) -> dict:
    pathway_dir = DATA_DIR / "pathways"
    pathway_dir.mkdir(parents=True, exist_ok=True)
    summary = []
    for pathway, frame in details.groupby("pathway_name", dropna=False):
        pathway = str(pathway)
        if not pathway or pathway == "nan":
            continue
        matrix = frame.pivot_table(
            index="source",
            columns="target",
            values="weighted_score",
            aggfunc="sum",
            fill_value=0.0,
        ).reindex(index=groups, columns=groups, fill_value=0.0)
        safe = "".join(char if char.isalnum() or char in "-_" else "_"
                       for char in pathway)
        matrix.to_csv(pathway_dir / f"{safe}_communication_matrix.csv")
        top = frame.nlargest(100, "weighted_score")
        top.to_csv(
            pathway_dir / f"{safe}_top_interactions.csv.gz",
            index=False,
            compression="gzip",
        )
        summary.append({
            "pathway": pathway,
            "n_interaction_records": int(len(frame)),
            "total_weighted_score": float(frame["weighted_score"].sum()),
            "top_sender": matrix.sum(axis=1).idxmax(),
            "top_receiver": matrix.sum(axis=0).idxmax(),
        })
    summary_df = pd.DataFrame(summary).sort_values(
        "total_weighted_score", ascending=False
    )
    summary_df.to_csv(DATA_DIR / "pathway_summary.csv", index=False)
    return {"n_pathways": int(len(summary_df))}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--permutations", type=int, default=N_PERMUTATIONS)
    parser.add_argument("--chunk-cells", type=int, default=512)
    parser.add_argument("--checkpoint-every", type=int, default=10)
    parser.add_argument("--force", action="store_true")
    parser.add_argument(
        "--skip-activity-export",
        action="store_true",
        help="Do not write the compact ligand/receptor activity parquet.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    ensure_directories()
    started = time.time()
    print("=" * 80)
    print("SCKNIFE ACTIVITY-BASED LIGAND-RECEPTOR COMMUNICATION")
    print("=" * 80)

    db, wanted_genes = load_interactions()
    parquet, node_names, cell_columns = activity_schema()
    metadata, groups = load_metadata(cell_columns)
    display_groups = groups
    print(
        f"Input: {parquet.metadata.num_rows:,} nodes x "
        f"{len(cell_columns):,} cells; {len(groups)} DE-derived final cell types"
    )

    activity, available_genes, activity_stats = stream_activity_subset(
        wanted_genes, node_names, cell_columns, args.chunk_cells
    )
    if not args.skip_activity_export:
        activity_df = pd.DataFrame(
            activity.T, index=available_genes, columns=cell_columns
        )
        activity_df.to_parquet(
            DATA_DIR / "ligand_receptor_activities.parquet",
            compression="zstd",
        )
        del activity_df
        gc.collect()

    group_to_code = {group: idx for idx, group in enumerate(groups)}
    labels = metadata[GROUP_COLUMN].map(group_to_code).to_numpy(dtype=np.int16)
    label_signature = hashlib.sha256(
        labels.tobytes() + json.dumps(groups, separators=(",", ":")).encode()
    ).hexdigest()
    observed_means, group_counts = group_means(activity, labels, len(groups))
    global_means = observed_means.mean(axis=1)
    specificity = observed_means / (global_means[:, None] + 1e-10)
    pairs = build_pairs(db, available_genes, observed_means)

    valid_db_rows = pairs["db_row"].nunique()
    print(
        f"Available nodes: {len(available_genes):,}; "
        f"valid DB rows: {valid_db_rows:,}; "
        f"LR components: {len(pairs):,}"
    )
    vectorization_error = validate_vectorization(
        observed_means, specificity, pairs
    )
    observed_array = communication_matrix(
        observed_means, specificity, pairs
    )
    observed = pd.DataFrame(
        observed_array, index=display_groups, columns=display_groups
    )

    details = interaction_details(
        observed_means, specificity, pairs, display_groups
    )
    details.to_csv(
        DATA_DIR / "interaction_details.csv.gz",
        index=False,
        compression="gzip",
    )
    pd.DataFrame(
        observed_means, index=available_genes, columns=display_groups
    ).to_csv(DATA_DIR / "mean_activity_by_celltype.csv.gz", compression="gzip")
    pd.DataFrame(
        specificity, index=available_genes, columns=display_groups
    ).to_csv(DATA_DIR / "specificity_by_celltype.csv.gz", compression="gzip")
    observed.to_csv(DATA_DIR / "communication_matrix.csv")

    pvalue_array, completed = run_permutations(
        activity=activity,
        labels=labels,
        observed=observed_array,
        specificity=specificity,
        pairs=pairs,
        n_groups=len(groups),
        n_permutations=args.permutations,
        checkpoint_every=args.checkpoint_every,
        force=args.force,
        label_signature=label_signature,
    )
    pvalues = pd.DataFrame(
        pvalue_array, index=display_groups, columns=display_groups
    )
    pvalues.to_csv(DATA_DIR / "communication_matrix_pvalues.csv")
    _, fdr_flat, _, _ = multipletests(
        pvalue_array.reshape(-1), method="fdr_bh"
    )
    fdr = pd.DataFrame(
        fdr_flat.reshape(pvalue_array.shape),
        index=display_groups,
        columns=display_groups,
    )
    fdr.to_csv(DATA_DIR / "communication_matrix_fdr.csv")
    observed.where(fdr < FDR_THRESHOLD, 0.0).to_csv(
        DATA_DIR / "communication_matrix_significant.csv"
    )

    centrality = network_centrality(observed)
    centrality.to_csv(DATA_DIR / "centrality_analysis.csv", index=False)

    meta_mask = metadata[METASTATIC_COLUMN].astype(str) == "E"
    meta_means, _ = group_means(activity[meta_mask.to_numpy()], labels[meta_mask], len(groups))
    nonmeta_means, _ = group_means(
        activity[(~meta_mask).to_numpy()], labels[~meta_mask], len(groups)
    )
    meta_comm = communication_matrix(meta_means, specificity, pairs)
    nonmeta_comm = communication_matrix(nonmeta_means, specificity, pairs)
    pd.DataFrame(
        meta_comm, index=display_groups, columns=display_groups
    ).to_csv(
        DATA_DIR / "communication_matrix_responder.csv"
    )
    pd.DataFrame(
        nonmeta_comm, index=display_groups, columns=display_groups
    ).to_csv(
        DATA_DIR / "communication_matrix_nonresponder.csv"
    )
    pd.DataFrame(
        meta_comm - nonmeta_comm,
        index=display_groups,
        columns=display_groups,
    ).to_csv(DATA_DIR / "dysregulated_signaling_responder_vs_nonresponder.csv")

    bx_mask = metadata[SITE_COLUMN].astype(str) == "Pre"
    surgery_mask = metadata[SITE_COLUMN].astype(str) == "On"
    bx_means, bx_counts = group_means(
        activity[bx_mask.to_numpy()], labels[bx_mask], len(groups)
    )
    surgery_means, surgery_counts = group_means(
        activity[surgery_mask.to_numpy()], labels[surgery_mask], len(groups)
    )
    bx_comm = communication_matrix(bx_means, specificity, pairs)
    surgery_comm = communication_matrix(surgery_means, specificity, pairs)
    pd.DataFrame(bx_comm, index=display_groups, columns=display_groups).to_csv(
        DATA_DIR / "communication_matrix_pre.csv"
    )
    pd.DataFrame(surgery_comm, index=display_groups, columns=display_groups).to_csv(
        DATA_DIR / "communication_matrix_on.csv"
    )
    pd.DataFrame(
        surgery_comm - bx_comm, index=display_groups, columns=display_groups
    ).to_csv(DATA_DIR / "dysregulated_signaling_on_vs_pre.csv")

    metadata[[GROUP_COLUMN, PATIENT_COLUMN, SITE_COLUMN, METASTATIC_COLUMN]].to_csv(DATA_DIR / "cell_metadata.csv.gz", compression="gzip")
    pd.DataFrame({
        GROUP_COLUMN: display_groups,
        "cell_count": group_counts,
    }).to_csv(DATA_DIR / "celltype_counts.csv", index=False)
    pathway_stats = save_pathway_outputs(details, display_groups)

    manifest = {
        "analysis": "Scknife activity-based ligand-receptor communication",
        "method_reference": "legacy step_11_cell_cell_signaling.py",
        "not_cellchat_probability_model": True,
        "inputs": {
            "activity_parquet": str(ACTIVITY_PATH),
            "final_celltype_h5ad": str(H5AD_PATH),
            "cellchat_db": str(CELLCHAT_DB_PATH),
        },
        "dimensions": {
            "activity_nodes": int(parquet.metadata.num_rows),
            "cells": len(cell_columns),
            "final_celltypes": len(display_groups),
            "available_signaling_nodes": len(available_genes),
            "valid_database_rows": int(valid_db_rows),
            "ligand_receptor_components": len(pairs),
        },
        "scknife_final_celltype_counts": {
            group: int(count)
            for group, count in zip(display_groups, group_counts)
        },
        "descriptive_contrasts": {
            "on_vs_pre": {
                "column": SITE_COLUMN,
                "direction": "On - Pre",
                "pre_cells": int(bx_mask.sum()),
                "on_cells": int(surgery_mask.sum()),
                "excluded_cells": 0,
                "pre_celltype_counts": dict(zip(display_groups, map(int, bx_counts))),
                "on_celltype_counts": dict(zip(display_groups, map(int, surgery_counts))),
                "inference": "descriptive cell-pooled; no patient-level test",
            }
        },
        "parameters": {
            "group_column": GROUP_COLUMN,
            "public_group_column": GROUP_COLUMN,
            "seed": RANDOM_SEED,
            "permutations": args.permutations,
            "completed_permutations": completed,
            "fdr_method": "Benjamini-Hochberg",
            "fdr_threshold": FDR_THRESHOLD,
            "confidence_weights": {"KEGG": 1.0, "other_evidence": 0.8, "missing": 0.5},
            "cofactor_multipliers": COFACTOR_MULTIPLIERS,
            "permutation_unit": "cell labels (legacy; not patient-blocked)",
            "specificity_during_permutation": "fixed observed specificity (legacy)",
            "label_signature": label_signature,
        },
        "activity_subset": activity_stats,
        "vectorized_vs_literal_max_abs_error": vectorization_error,
        "pathways": pathway_stats,
        "software": software_versions(),
        "elapsed_seconds": time.time() - started,
    }
    write_json(DATA_DIR / "run_manifest.json", manifest)
    (DATA_DIR / "RUN_COMPLETE").touch()
    print(f"Complete in {manifest['elapsed_seconds'] / 60:.1f} minutes")


if __name__ == "__main__":
    main()
