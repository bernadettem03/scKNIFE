import os
import pandas as pd
import networkx as nx
from scknife.utils.helpers import add_edge_with_source


def load_cell_state(csv_path):
    print("  [cell_state] Loading cellMarker 2.0 (cancer breast tissue)...")
    G = nx.Graph()
    if os.path.exists(csv_path):
        try:
            cm2_csv = pd.read_csv(csv_path)
            if "Cell name" in cm2_csv.columns and "Cell marker" in cm2_csv.columns:
                cm2_df = cm2_csv[["Cell name", "Cell marker"]].copy()
                cm2_df.columns = ["cell_state", "genesymbol"]
                cm2_df = cm2_df.dropna()
                for _, row in cm2_df.iterrows():
                    g = str(row["genesymbol"])
                    s = str(row["cell_state"])
                    if not G.has_node(g):
                        G.add_node(g, kind="gene")
                    if not G.has_node(s):
                        G.add_node(s, kind="cell_state")
                    add_edge_with_source(G, g, s, source="cellmarker2.0", type="gene-cellstate")
                n_genes = sum(1 for _, a in G.nodes(data=True) if a.get("kind") == "gene")
                n_states = sum(1 for _, a in G.nodes(data=True) if a.get("kind") == "cell_state")
                print(f"    cellMarker2.0 loaded: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
                print(f"      (genes={n_genes}, cell_states={n_states})")
                return {"graph": G, "node_kinds": {"gene": n_genes, "cell_state": n_states}}
            else:
                print("    cellMarker2.0 CSV columns not recognized")
                return {"graph": nx.Graph(), "node_kinds": {}}
        except Exception as e:
            print(f"    Could not load cellMarker2.0 CSV: {str(e)[:60]}")
            return {"graph": nx.Graph(), "node_kinds": {}}
    else:
        print(f"    cellMarker2.0 CSV not found: {csv_path}")
        return {"graph": nx.Graph(), "node_kinds": {}}
