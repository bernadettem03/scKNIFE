import os
import pandas as pd
import networkx as nx
from scknife.utils.helpers import add_edge_with_source


def load_cellchat(csv_path):
    print("  [cellchat] Loading...")
    if not os.path.exists(csv_path):
        print("  XX File not found:", csv_path)
        return {"graph": nx.Graph(), "node_kinds": {}}
    interaction = pd.read_csv(csv_path)
    keywords = ["Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact"]
    filtered = interaction[interaction["annotation"].isin(keywords)]
    G = nx.DiGraph()
    for _, row in filtered.iterrows():
        ligand = row["ligand.symbol"]
        if not G.has_node(ligand):
            G.add_node(ligand, kind="cellchat_entity")
        for rec in str(row["receptor.symbol"]).split(","):
            rec = rec.strip()
            if rec:
                if not G.has_node(rec):
                    G.add_node(rec, kind="cellchat_entity")
                add_edge_with_source(G, ligand, rec, source="cellchat", pathway=row.get("pathway_name"))
    print(f"    {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
    return {"graph": G, "node_kinds": {"cellchat_entity": G.number_of_nodes()}}
