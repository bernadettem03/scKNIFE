#!/usr/bin/env python3
"""Validate the frozen scKNIFE graph releases against their audit manifest."""

from __future__ import annotations

import hashlib
import json
from collections import Counter
from pathlib import Path

import networkx as nx


ROOT = Path(__file__).resolve().parents[1]
MANIFEST_PATH = ROOT / "graphs" / "manifest.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def source_contains(value: object, source: str) -> bool:
    return source.lower() in str(value or "").lower()


def biological_program_counts(graph: nx.Graph) -> dict[str, int]:
    kinds = Counter(attrs.get("kind", "?") for _, attrs in graph.nodes(data=True))
    reaction_programs = {
        node
        for left, right, attrs in graph.edges(data=True)
        if source_contains(attrs.get("source"), "task_rxn")
        for node in (left, right)
        if graph.nodes[node].get("kind") == "reaction"
    }
    cellchat_pathways = {
        attrs["pathway"]
        for _, _, attrs in graph.edges(data=True)
        if source_contains(attrs.get("source"), "cellchat") and attrs.get("pathway")
    }
    counts = {
        "hallmark_gene_sets": kinds["geneset"],
        "metabolic_tasks": kinds["task"],
        "reaction_associated_programs": len(reaction_programs),
        "cell_type_marker_sets": kinds["cell_type"],
        "cell_state_programs": kinds["cell_state"],
        "go_biological_processes": kinds["GO_term"],
        "human_gem_pathways": kinds["pathway"],
        "cellchat_pathways": len(cellchat_pathways),
    }
    counts["total"] = sum(counts.values())
    return counts


def canonical_edges(graph: nx.Graph) -> set[frozenset[str]]:
    return {frozenset((str(left), str(right))) for left, right in graph.edges()}


def validate_version(name: str, expected: dict[str, object]) -> nx.Graph:
    path = ROOT / str(expected["file"])
    assert path.is_file(), f"Missing {name} graph: {path}"
    assert path.stat().st_size == expected["bytes"], f"Unexpected {name} file size"
    assert sha256(path) == expected["sha256"], f"Unexpected {name} SHA-256"

    graph = nx.read_graphml(path)
    assert graph.number_of_nodes() == expected["nodes"], f"Unexpected {name} node count"
    assert graph.number_of_edges() == expected["edges"], f"Unexpected {name} edge count"
    assert nx.number_connected_components(graph) == expected["connected_components"], f"Unexpected {name} component count"
    kinds = Counter(attrs.get("kind", "?") for _, attrs in graph.nodes(data=True))
    assert dict(kinds) == expected["node_kinds"], f"Unexpected {name} node-kind counts"
    assert biological_program_counts(graph) == expected["biological_programs"], f"Unexpected {name} biological-program counts"
    print(f"PASS {name}: {graph.number_of_nodes():,} nodes, {graph.number_of_edges():,} edges, {expected['biological_programs']['total']:,} programs")
    return graph


def main() -> None:
    manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    expected = manifest["versions"]
    v1 = validate_version("version_1", expected["version_1"])
    v2 = validate_version("version_2", expected["version_2"])

    assert set(v1.nodes()).issubset(v2.nodes()), "Version 1 nodes are not preserved in Version 2"
    assert canonical_edges(v1).issubset(canonical_edges(v2)), "Version 1 edges are not preserved in Version 2"
    increment = expected["version_2"]["increment_over_version_1"]
    assert v2.number_of_nodes() - v1.number_of_nodes() == increment["nodes"]
    assert v2.number_of_edges() - v1.number_of_edges() == increment["edges"]
    print("PASS version_2 preserves every version_1 node and edge")


if __name__ == "__main__":
    main()
