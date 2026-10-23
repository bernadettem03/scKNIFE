import os
import pandas as pd
import networkx as nx
from scknife.utils.helpers import add_edge_with_source


def load_breast_ppi(ppi_path, tpm_threshold=5, intconf_threshold=0.9):
    print("  [breast_ppi] Loading...")
    if not os.path.exists(ppi_path):
        print("  XX File not found:", ppi_path)
        return {"graph": nx.Graph(), "node_kinds": {}}

    ppi = pd.read_csv(ppi_path, sep="\t", header=None,
                      names=["gene1", "gene2", "intconf", "tpm1", "tpm2"])
    ppi_f = ppi[(ppi["intconf"] > intconf_threshold) & (ppi["tpm1"] >= tpm_threshold) & (ppi["tpm2"] >= tpm_threshold)]
    print(f"    {len(ppi_f)} interactions (filtered from {len(ppi)})")

    G = nx.Graph()
    for _, r in ppi_f.iterrows():
        g1, g2 = r["gene1"], r["gene2"]
        if not G.has_node(g1):
            G.add_node(g1, kind="gene", source="breast_PPI")
        if not G.has_node(g2):
            G.add_node(g2, kind="gene", source="breast_PPI")
        add_edge_with_source(G, g1, g2, source="breast_PPI", intconf=r["intconf"],
                             tpm_gene1=r["tpm1"], tpm_gene2=r["tpm2"])

    print(f"    {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
    return {"graph": G, "node_kinds": {"gene": G.number_of_nodes()}}
