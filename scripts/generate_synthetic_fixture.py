"""Generate invented counts for software regression tests; no patient data are used."""
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp
import networkx as nx

ROOT = Path(__file__).resolve().parents[1]

def generate():
    rng = np.random.default_rng(42)
    graph = nx.read_graphml(ROOT / "graphs/version_1/unified_graph.graphml")
    genes = sorted(n for n, a in graph.nodes(data=True) if a.get("kind") == "gene")
    chosen = sorted(rng.choice(genes, 1000, replace=False))
    groups = np.repeat(np.arange(4), 250)
    rng.shuffle(groups)
    profiles = rng.gamma(1.5, 0.8, (4, 1000))
    library = rng.lognormal(0, 0.3, (1000, 1))
    counts = rng.poisson(profiles[groups] * library).astype(np.int32)
    obs = pd.DataFrame({"synthetic_group": [f"Synthetic group {g+1}" for g in groups]},
                       index=[f"synthetic_cell_{i:04d}" for i in range(1000)])
    data = ad.AnnData(sp.csr_matrix(counts), obs=obs, var=pd.DataFrame(index=chosen))
    data.uns["data_origin"] = "Fully synthetic counts and group labels; no patient data"
    data.uns["generation_seed"] = 42
    out = ROOT / "tutorial/data"
    data.write_h5ad(out / "synthetic_regression_1000.h5ad", compression="gzip")
    obs.to_csv(out / "synthetic_regression_1000_metadata.csv")
    return data

if __name__ == "__main__":
    generate()
