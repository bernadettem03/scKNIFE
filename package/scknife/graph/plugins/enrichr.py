import os
import glob
import pandas as pd
import networkx as nx
from scknife.utils.helpers import add_edge_with_source


def _extract_enrichr_basename(fpath):
    bn = os.path.basename(fpath)
    bn = bn.removesuffix(".csv").removesuffix(".CSV")
    if " Edges" in bn:
        bn = bn.split(" Edges")[0]
    if " Gene Edges" in bn:
        bn = bn.split(" Gene Edges")[0]
    if bn.endswith(" Gene"):
        bn = bn.removesuffix(" Gene")
    if "." in bn and ("GO_Biological" in bn or "Reactome" in bn or "TRRUST" in bn or "Tabula_Muris" in bn):
        bn = bn.split(".")[0]
    return bn


def _load_enrichr_file(fpath, kind_map, default_kind="knowledge_base"):
    bn = _extract_enrichr_basename(fpath)
    edge_kind = kind_map.get(bn, default_kind)
    is_digraph = "TF" in bn

    edges = pd.read_csv(fpath)[["source_label", "target_label"]].dropna()
    edges = edges.rename(columns={"source_label": "source", "target_label": "target"})
    G = nx.DiGraph() if is_digraph else nx.Graph()

    for _, row in edges.iterrows():
        src_n, tgt_n = row["source"], row["target"]
        if not G.has_node(src_n):
            G.add_node(src_n, kind=edge_kind)
        elif "kind" not in G.nodes[src_n]:
            G.nodes[src_n]["kind"] = edge_kind
        if not G.has_node(tgt_n):
            G.add_node(tgt_n, kind="gene")
        elif "kind" not in G.nodes[tgt_n]:
            G.nodes[tgt_n]["kind"] = "gene"
        add_edge_with_source(G, src_n, tgt_n, source=bn)
    return G


def _find_enrichr_file(enrich_dir, keyword):
    matches = [f for f in glob.glob(os.path.join(enrich_dir, "*"))
               if (f.endswith("Edges.csv") or f.endswith("edges.csv"))
               and keyword in f and "Human Gene Edges" not in f]
    return matches[0] if matches else None


ENRICHR_KIND_MAP = {
    "GO_Biological_Process_2021":          "GO_term",
    "Reactome_2022":                       "Reactome_pathway",
    "TRRUST_Transcription_Factors_2019":   "TF",
}


def _load_enrichr_plugin(enrich_dir, keyword, plugin_name, kind_counts):
    print(f"  [enrichr_{plugin_name}] Loading...")
    fpath = _find_enrichr_file(enrich_dir, keyword)
    if not fpath:
        return {"graph": nx.Graph(), "node_kinds": {}}
    G = _load_enrichr_file(fpath, ENRICHR_KIND_MAP)
    n_src = sum(1 for _, a in G.nodes(data=True) if a.get("kind") == kind_counts["src_kind"])
    n_gn = sum(1 for _, a in G.nodes(data=True) if a.get("kind") == "gene")
    print(f"    {G.number_of_nodes()} nodes, {G.number_of_edges()} edges ({kind_counts['src_label']}={n_src}, genes={n_gn})")
    return {"graph": G, "node_kinds": {kind_counts["src_kind"]: n_src, "gene": n_gn}}


def load_enrichr_go(enrich_dir):
    return _load_enrichr_plugin(enrich_dir, "GO_Biological", "go",
                                {"src_kind": "GO_term", "src_label": "GO_terms"})


def load_enrichr_reactome(enrich_dir):
    return _load_enrichr_plugin(enrich_dir, "Reactome", "reactome",
                                {"src_kind": "Reactome_pathway", "src_label": "pathways"})


def load_enrichr_kegg(enrich_dir):
    return _load_enrichr_plugin(enrich_dir, "KEGG", "kegg",
                                {"src_kind": "KEGG_pathway", "src_label": "pathways"})


def load_enrichr_trrust(enrich_dir):
    return _load_enrichr_plugin(enrich_dir, "TRRUST", "trrust",
                                {"src_kind": "TF", "src_label": "TFs"})


def load_enrichr_archs4(enrich_dir):
    return _load_enrichr_plugin(enrich_dir, "ARCHS4", "archs4",
                                {"src_kind": "TF", "src_label": "TFs"})


def load_enrichr_tabula(enrich_dir):
    return _load_enrichr_plugin(enrich_dir, "Tabula_Muris", "tabula",
                                {"src_kind": "cell_type", "src_label": "cell_types"})
