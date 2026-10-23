import os
import ast
import pandas as pd
import networkx as nx
from scknife.utils.helpers import add_edge_with_source


def load_hallmark(csv_path):
    print("  [hallmark] Loading...")
    G = nx.Graph()
    try:
        msig = pd.read_csv(csv_path)
        msig_h = msig[msig["standard_name"].str.upper().str.contains("HALLMARK")].copy()
        msig_h["genes"] = msig_h["GENE_SYMBOLS"].apply(
            lambda s: ast.literal_eval(s) if isinstance(s, str) and s.startswith("[") else []
        )
        for _, row in msig_h.iterrows():
            sn = row["standard_name"]
            G.add_node(sn, bipartite=0, kind="geneset")
            for g in row["genes"]:
                G.add_node(g, bipartite=1, kind="gene")
                add_edge_with_source(G, sn, g, source="hallmark")
    except Exception as e:
        print("  XX", e)
    n_gs = sum(1 for _, a in G.nodes(data=True) if a.get("kind") == "geneset")
    n_g = sum(1 for _, a in G.nodes(data=True) if a.get("kind") == "gene")
    print(f"    {G.number_of_nodes()} nodes, {G.number_of_edges()} edges (genesets={n_gs}, genes={n_g})")
    return {"graph": G, "node_kinds": {"geneset": n_gs, "gene": n_g}}
