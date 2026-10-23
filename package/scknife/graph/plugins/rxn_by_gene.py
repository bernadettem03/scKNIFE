import os
import pandas as pd
import networkx as nx
from scknife.utils.helpers import add_edge_with_source


def load_rxn_by_gene(csv_path):
    print("  [rxn_by_gene] Loading...")
    if not os.path.exists(csv_path):
        print("  XX File not found:", csv_path)
        return {"graph": nx.Graph(), "node_kinds": {}}
    df = pd.read_csv(csv_path, index_col=0)
    G = nx.Graph()
    G.add_nodes_from(df.index, bipartite="reaction", kind="reaction")
    G.add_nodes_from(df.columns, bipartite="gene", kind="gene")
    for rxn in df.index:
        for gene in df.columns:
            if df.loc[rxn, gene] != 0:
                add_edge_with_source(G, rxn, gene, source="rxn_gene")
    print(f"    {G.number_of_nodes()} nodes, {G.number_of_edges()} edges (reactions={len(df.index)}, genes={len(df.columns)})")
    return {"graph": G, "node_kinds": {"reaction": len(df.index), "gene": len(df.columns)}}
