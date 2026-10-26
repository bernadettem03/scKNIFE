import os
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
import networkx as nx
from scknife.utils.helpers import add_edge_with_source, sanitize_name


SBML_NS = "http://www.sbml.org/sbml/level3/version1/core"
FBC_NS = "http://www.sbml.org/sbml/level3/version1/fbc/version2"
GRP_NS = "http://www.sbml.org/sbml/level3/version1/groups/version1"
COMPARTMENT_NAMES = {
    "c": "cytosol", "m": "mitochondria", "e": "extracellular",
    "n": "nucleus", "g": "golgi", "l": "lysosome",
    "i": "inner_mitochondria", "x": "peroxisome", "r": "endoplasmic_reticulum",
}


def load_human_gem(xml_path):
    print("  [human_gem] Parsing...")
    if not os.path.exists(xml_path):
        print("  XX File not found:", xml_path)
        return {"graph": nx.Graph(), "node_kinds": {}}

    tree = ET.parse(xml_path)
    root = tree.getroot()
    model = root.find(f"{{{SBML_NS}}}model")

    G = nx.Graph()
    stats = Counter()
    met_id_to_name = {}
    rxn_id_to_name = {}

    met_name_to_compartments = defaultdict(set)
    sp_el = model.find(f"{{{SBML_NS}}}listOfSpecies")
    if sp_el is not None:
        for sp in sp_el.findall(f"{{{SBML_NS}}}species"):
            sid = sp.get("id")
            sname = sp.get("name", sid)
            comp = sp.get("compartment", "")
            met_name_to_compartments[sname].add(comp)
            readable_name = sanitize_name(sname)
            met_id_to_name[sid] = readable_name
            if not G.has_node(readable_name):
                G.add_node(readable_name, kind="metabolite",
                           original_id=sid, name=sname, source="Human_GEM")
            stats["metabolites"] += 1
    for met_name, comps in met_name_to_compartments.items():
        rname = sanitize_name(met_name)
        if G.has_node(rname):
            G.nodes[rname]["compartments"] = {COMPARTMENT_NAMES.get(c, c) for c in comps}

    ensembl_to_symbol = {}
    gp_list = model.find(f"{{{FBC_NS}}}listOfGeneProducts")
    if gp_list is not None:
        for gp in gp_list.findall(f"{{{FBC_NS}}}geneProduct"):
            gid = gp.get(f"{{{FBC_NS}}}id")
            label = gp.get(f"{{{FBC_NS}}}label", gid)
            ensembl_to_symbol[gid] = label
            if not G.has_node(label):
                G.add_node(label, kind="gene", gene_symbol=label,
                           original_id=gid, ensembl_id=gid, source="Human_GEM")
            stats["genes"] += 1

    rxn_el = model.find(f"{{{SBML_NS}}}listOfReactions")
    if rxn_el is not None:
        for rxn in rxn_el.findall(f"{{{SBML_NS}}}reaction"):
            rid = rxn.get("id")
            rname = rxn.get("name", rid)
            rev = rxn.get("reversible", "false")
            readable_name = sanitize_name(rname)
            rxn_id_to_name[rid] = readable_name
            G.add_node(readable_name, kind="reaction",
                       original_id=rid, name=rname, reversible=rev, source="Human_GEM")
            stats["reactions"] += 1

            rel = rxn.find(f"{{{SBML_NS}}}listOfReactants")
            if rel is not None:
                for sr in rel.findall(f"{{{SBML_NS}}}speciesReference"):
                    met_sid = sr.get("species")
                    sto = float(sr.get("stoichiometry", "1"))
                    if met_sid:
                        met_readable = met_id_to_name.get(met_sid, met_sid)
                        add_edge_with_source(G, readable_name, met_readable,
                            source="Human_GEM", edge_type="reaction-metabolite",
                            role="reactant", stoichiometry=sto)
                        stats["rxn-met edges"] += 1
            pel = rxn.find(f"{{{SBML_NS}}}listOfProducts")
            if pel is not None:
                for sr in pel.findall(f"{{{SBML_NS}}}speciesReference"):
                    met_sid = sr.get("species")
                    sto = float(sr.get("stoichiometry", "1"))
                    if met_sid:
                        met_readable = met_id_to_name.get(met_sid, met_sid)
                        add_edge_with_source(G, readable_name, met_readable,
                            source="Human_GEM", edge_type="reaction-metabolite",
                            role="product", stoichiometry=sto)
                        stats["rxn-met edges"] += 1
            gpa = rxn.find(f"{{{FBC_NS}}}geneProductAssociation")
            if gpa is not None:
                for ref in gpa.findall(f".//{{{FBC_NS}}}geneProductRef"):
                    gid = ref.get(f"{{{FBC_NS}}}geneProduct")
                    if gid:
                        gsym = ensembl_to_symbol.get(gid, gid)
                        add_edge_with_source(G, readable_name, gsym,
                            source="Human_GEM", edge_type="gene-reaction", ensembl_id=gid)
                        stats["gene-rxn edges"] += 1

    grp_list = model.find(f"{{{GRP_NS}}}listOfGroups")
    if grp_list is not None:
        for grp in grp_list.findall(f"{{{GRP_NS}}}group"):
            gid = grp.get(f"{{{GRP_NS}}}id", "")
            gname = grp.get(f"{{{GRP_NS}}}name", "").strip()
            if not gid and not gname:
                continue
            readable_name = sanitize_name(gname) if gname else sanitize_name(gid)
            G.add_node(readable_name, kind="pathway",
                       original_id=gid, name=gname, source="Human_GEM")
            stats["pathways"] += 1
            for m in grp.findall(f".//{{{GRP_NS}}}member"):
                rid = m.get(f"{{{GRP_NS}}}idRef", "")
                if rid:
                    rxn_readable = rxn_id_to_name.get(rid, rid)
                    add_edge_with_source(G, readable_name, rxn_readable,
                        source="Human_GEM", edge_type="pathway-reaction")
                    stats["pathway-rxn edges"] += 1

    for gsym in list(ensembl_to_symbol.values()):
        if not G.has_node(gsym):
            continue
        rxn_neighbors = [n for n in G.neighbors(gsym) if G.nodes[n].get("kind") == "reaction"]
        for rxn in rxn_neighbors:
            met_neighbors = [n for n in G.neighbors(rxn) if G.nodes[n].get("kind") == "metabolite"]
            for met in met_neighbors:
                if not G.has_edge(gsym, met):
                    add_edge_with_source(G, gsym, met,
                        source="Human_GEM", edge_type="gene-metabolite", via_reaction=rxn)
                    stats["gene-met shortcut edges"] += 1

    print(f"    mets={stats['metabolites']} genes={stats['genes']} rxns={stats['reactions']} paths={stats['pathways']} edges={G.number_of_edges()}")

    n_gene = sum(1 for _, a in G.nodes(data=True) if a.get("kind") == "gene")
    n_met = sum(1 for _, a in G.nodes(data=True) if a.get("kind") == "metabolite")
    n_rxn = sum(1 for _, a in G.nodes(data=True) if a.get("kind") == "reaction")
    n_pw = sum(1 for _, a in G.nodes(data=True) if a.get("kind") == "pathway")

    return {"graph": G, "node_kinds": {"gene": n_gene, "metabolite": n_met, "reaction": n_rxn, "pathway": n_pw}}
