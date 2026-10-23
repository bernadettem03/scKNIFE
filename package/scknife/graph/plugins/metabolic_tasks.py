import os
import pandas as pd
import networkx as nx
from scknife.utils.helpers import add_edge_with_source


def load_metabolic_tasks(csv_path):
    print("  [metabolic_tasks] Loading...")
    if not os.path.exists(csv_path):
        print("  XX File not found:", csv_path)
        return {"graph": nx.Graph(), "node_kinds": {}}
    df = pd.read_csv(csv_path, index_col=0)
    G = nx.Graph()
    G.add_nodes_from(df.index, bipartite="task", kind="task")
    G.add_nodes_from(df.columns, bipartite="gene", kind="gene")
    for rxn in df.index:
        for gene in df.columns:
            if df.loc[rxn, gene] != 0:
                add_edge_with_source(G, rxn, gene, source="metabolic_task")
    print(f"    {G.number_of_nodes()} nodes, {G.number_of_edges()} edges (tasks={len(df.index)}, genes={len(df.columns)})")
    return {"graph": G, "node_kinds": {"task": len(df.index), "gene": len(df.columns)}}
