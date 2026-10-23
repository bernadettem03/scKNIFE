import numpy as np
import scipy.sparse as sp


def build_Rg_mapping(adata, G):
    print("Building gene-to-KG node mapping (Rg)...")
    gene_names = adata.var_names.tolist()
    node_list = list(G.nodes)
    node_name_to_index = {n: i for i, n in enumerate(node_list)}
    mapping_idx = []

    for i, g in enumerate(gene_names):
        j = node_name_to_index.get(g, None)
        if j is not None:
            mapping_idx.append((i, j))

    rows, cols = zip(*mapping_idx) if mapping_idx else ([], [])
    Rg = sp.csr_matrix((np.ones(len(mapping_idx)), (rows, cols)),
                       shape=(len(gene_names), len(node_list)))
    adata.uns['Rg'] = Rg
    adata.uns['node_names'] = node_list
    print(f"Mapped {len(mapping_idx)} genes to KG nodes.")
    return Rg, mapping_idx, node_list
