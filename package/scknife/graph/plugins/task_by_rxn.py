import os
import pandas as pd
import networkx as nx
from scknife.utils.helpers import add_edge_with_source


def load_task_by_rxn(csv_path):
    print("  [task_by_rxn] Loading...")
    if not os.path.exists(csv_path):
        print("  XX File not found:", csv_path)
        return {"graph": nx.Graph(), "node_kinds": {}}
    df = pd.read_csv(csv_path, index_col=0)
    G = nx.Graph()
    G.add_nodes_from(df.index, bipartite="task", kind="task")
    G.add_nodes_from(df.columns, bipartite="reaction", kind="reaction")
    for task in df.index:
        for rxn in df.columns:
            if df.loc[task, rxn] != 0:
                add_edge_with_source(G, task, rxn, source="task_rxn")
    print(f"    {G.number_of_nodes()} nodes, {G.number_of_edges()} edges (tasks={len(df.index)}, reactions={len(df.columns)})")
    return {"graph": G, "node_kinds": {"task": len(df.index), "reaction": len(df.columns)}}
