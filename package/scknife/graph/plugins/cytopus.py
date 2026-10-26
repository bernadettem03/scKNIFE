import networkx as nx
from scknife.utils.helpers import add_edge_with_source


def _dict_to_graph(mapping, source_kind):
    G = nx.Graph()
    for src, genes in mapping.items():
        if not src:
            continue
        G.add_node(src, kind=source_kind)
        for g in genes:
            if not g:
                continue
            if not G.has_node(g):
                G.add_node(g, kind="gene")
            add_edge_with_source(G, src, g, source=f"cytopus_{source_kind}")
    return G


def load_cytopus_processes():
    import cytopus as cp
    print("  [cytopus_processes] Loading...")
    A = cp.KnowledgeBase()
    G = _dict_to_graph(A.processes, "process")
    n_proc = sum(1 for _, a in G.nodes(data=True) if a.get("kind") == "process")
    n_gene = sum(1 for _, a in G.nodes(data=True) if a.get("kind") == "gene")
    print(f"    {G.number_of_nodes()} nodes, {G.number_of_edges()} edges (processes={n_proc}, genes={n_gene})")
    return {"graph": G, "node_kinds": {"process": n_proc, "gene": n_gene}}


def load_cytopus_identities():
    import cytopus as cp
    print("  [cytopus_identities] Loading...")
    A = cp.KnowledgeBase()
    G = _dict_to_graph(A.identities, "identity")
    n_id = sum(1 for _, a in G.nodes(data=True) if a.get("kind") == "identity")
    n_gene = sum(1 for _, a in G.nodes(data=True) if a.get("kind") == "gene")
    print(f"    {G.number_of_nodes()} nodes, {G.number_of_edges()} edges (identities={n_id}, genes={n_gene})")
    return {"graph": G, "node_kinds": {"identity": n_id, "gene": n_gene}}
