import os
import networkx as nx
from collections import Counter

from scknife.utils.helpers import add_edge_with_source


def merge_graph(G_unified, G_source, source_name):
    added_nodes = 0
    added_edges = 0
    for n, attrs in G_source.nodes(data=True):
        if not G_unified.has_node(n):
            G_unified.add_node(n, **attrs)
            added_nodes += 1
    for u, v, attrs in G_source.edges(data=True):
        if not G_unified.has_edge(u, v):
            G_unified.add_edge(u, v, **attrs)
            added_edges += 1
    print(f"   + {source_name:30s}  +{added_nodes:>7,} nodes  +{added_edges:>7,} edges")
    return added_nodes, added_edges


def unify_plugins(plugins, graph_name="unified"):
    G = nx.Graph()
    total_n, total_e = 0, 0
    for pname, result in plugins.items():
        graph_obj = result if isinstance(result, nx.Graph) else result.get("graph", nx.Graph())
        if graph_obj.number_of_nodes() == 0:
            print(f"   - {pname:30s}  SKIP (empty)")
            continue
        an, ae = merge_graph(G, graph_obj, pname)
        total_n += an
        total_e += ae
    print(f"   {'-' * 50}")
    print(f"   TOTAL: {G.number_of_nodes():,} nodes  {G.number_of_edges():,} edges")
    return G


def plugin_report(results, title="GRAPH PLUGINS"):
    print()
    print("=" * 70)
    print(f"  {title}")
    print("=" * 70)
    print(f"  {'Plugin':30s} {'Nodes':>8s} {'Edges':>8s}  {'Kinds'}")
    print(f"  {'-' * 65}")
    total_n, total_e = 0, 0
    for name, result in results.items():
        graph_obj = result if isinstance(result, nx.Graph) else result.get("graph", nx.Graph())
        n = graph_obj.number_of_nodes()
        e = graph_obj.number_of_edges()
        nk = result.get("node_kinds", {}) if isinstance(result, dict) else {}
        kinds = ", ".join(f"{k}={v:,}" for k, v in sorted(nk.items(), key=lambda x: -x[1])) if nk else ""
        print(f"  {name:30s} {n:>8,} {e:>8,}  {kinds}")
        total_n += n
        total_e += e
    print(f"  {'-' * 65}")
    print(f"  {'TOTAL (before dedup)':30s} {total_n:>8,} {total_e:>8,}")


def construct_unified_graph(plugins, graph_name="unified"):
    return unify_plugins(plugins, graph_name)
