import re
import os
import time
import json
from contextlib import contextmanager
import numpy as np
import pandas as pd
import networkx as nx


def get_tuning_cells(n_cells, n_graph_nodes):
    n_tune = int(0.1 * np.sqrt(n_cells * n_graph_nodes))
    return min(20000, max(5000, n_tune))


def _fmt_duration(seconds):
    h, r = divmod(seconds, 3600)
    m, s = divmod(r, 60)
    if h >= 1: return f"{int(h)}h {int(m)}m {s:.0f}s"
    if m >= 1: return f"{int(m)}m {s:.0f}s"
    return f"{s:.1f}s"


def format_duration(seconds):
    h, r = divmod(seconds, 3600)
    m, s = divmod(r, 60)
    if h >= 1: return f"{int(h)}h {int(m)}m {s:.2f}s"
    if m >= 1: return f"{int(m)}m {s:.2f}s"
    return f"{s:.2f}s"


class ProgressTracker:
    def __init__(self, total):
        self.total = total
        self.current = 0
        self.t0 = None

    def start(self):
        self.t0 = time.time()

    def step_done(self, label=""):
        self.current += 1
        elapsed = time.time() - self.t0
        pct = self.current / self.total * 100
        etc = (elapsed / self.current) * (self.total - self.current)
        bar_len = 25
        filled = int(bar_len * self.current / self.total)
        bar = "\u2588" * filled + "\u2591" * (bar_len - filled)
        print(f"\n  Progress: [{bar}] {pct:.0f}%  |  Elapsed: {_fmt_duration(elapsed)}  |  ETC: {_fmt_duration(etc)}  |  {label}\n")


class TeeOutput:
    def __init__(self, *files):
        self.files = files
    def write(self, data):
        for f in self.files:
            f.write(data)
            f.flush()
    def flush(self):
        for f in self.files:
            f.flush()


def sanity_summary(step_name, checks):
    for label, passed, detail in checks:
        icon = "\u2713" if passed else "\u2717"
        print(f"    {icon} {label}: {detail}")


@contextmanager
def timed_step(label, step_timings, progress=None, step_name=""):
    t0 = time.time()
    print(label)
    try:
        yield
    finally:
        elapsed = time.time() - t0
        step_timings[label] = elapsed
        print(f"  \u23f1  Step runtime: {format_duration(elapsed)} ({elapsed:.2f}s)")
        if progress is not None:
            progress.step_done(step_name)


def store_H_to_obs(adata, H, prefix):
    print(f"Storing {prefix} components in adata.obs...")
    H_df = pd.DataFrame(H.T, index=adata.obs_names, columns=[f"{prefix}_{i}" for i in range(H.shape[0])])
    old = [c for c in adata.obs.columns if str(c).startswith(prefix + "_")]
    adata.obs = pd.concat([adata.obs.drop(columns=old), H_df], axis=1, copy=False)


def h_projection(W_graph, Rg, H_train, X_test, train_idx, test_idx,
                 adata_orig, k, batch_size=10000):
    import gc
    import scipy.sparse as sp
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        from scipy.optimize import nnls, lsq_linear
    import anndata

    n_test = X_test.shape[1]
    n_all = adata_orig.n_obs

    # Convert to float32 to halve memory
    if X_test.dtype != np.float32:
        X_test = X_test.astype(np.float32)
    if H_train.dtype != np.float32:
        H_train = H_train.astype(np.float32)

    basis = (Rg @ W_graph)
    if sp.issparse(basis):
        basis = basis.toarray().astype(np.float32)
    else:
        basis = basis.astype(np.float32)

    print(f"NNLS projecting {n_test:,} test cells onto trained basis ({k} components)...")
    H_test = np.zeros((k, n_test), dtype=np.float32)
    for i in range(n_test):
        try:
            h, _ = nnls(basis, X_test[:, i], maxiter=300000)
        except RuntimeError:
            try:
                res = lsq_linear(basis, X_test[:, i], bounds=(0, np.inf),
                                 method='trf', max_iter=1000, tol=1e-6)
                h = res.x
            except Exception:
                h = np.zeros(k, dtype=np.float32)
        H_test[:, i] = h.astype(np.float32)
        if (i + 1) % 5000 == 0:
            print(f"    Projected {i+1:,}/{n_test:,} cells")

    del X_test
    gc.collect()

    H_combined = np.zeros((k, n_all), dtype=np.float32)
    H_combined[:, train_idx] = H_train
    H_combined[:, test_idx] = H_test
    del H_train, H_test
    gc.collect()

    # Reconstruct expression in chunks to limit peak memory
    n_genes = basis.shape[0]
    X_recon = np.empty((n_all, n_genes), dtype=np.float32)
    for start in range(0, n_all, batch_size):
        end = min(start + batch_size, n_all)
        X_recon[start:end] = (basis @ H_combined[:, start:end]).T
        print(f"    Recon batch {start//batch_size + 1}/{(n_all + batch_size - 1)//batch_size}: cells {start}-{end}")
    del basis
    gc.collect()

    adata_out = anndata.AnnData(
        X=X_recon,
        obs=adata_orig.obs.copy(),
        var=adata_orig.var.copy()
    )
    del X_recon
    gc.collect()

    store_H_to_obs(adata_out, H_combined, "H_graph")
    adata_out.uns['W_graph'] = W_graph
    adata_out.uns['Rg'] = Rg
    adata_out.uns['node_names'] = adata_orig.uns.get('node_names', [])
    adata_out.uns['H_combined'] = H_combined
    print(f"  Assembled: {adata_out.n_obs:,} cells, {adata_out.n_vars:,} genes")
    return adata_out


def add_edge_with_source(G, u, v, source, **attrs):
    if G.has_edge(u, v):
        existing = G[u][v].get("source")
        if isinstance(existing, set):
            existing.add(source)
        else:
            G[u][v]["source"] = {existing, source} if existing else {source}
        for k, val in attrs.items():
            G[u][v][k] = val
    else:
        attrs["source"] = source
        G.add_edge(u, v, **attrs)


def sanitize_name(name):
    name = str(name)
    for ch in " :/()[]-+.,;'\"\"":
        name = name.replace(ch, "_")
    while "__" in name:
        name = name.replace("__", "_")
    return name.strip("_")


def make_graphml_safe(obj):
    if isinstance(obj, (str, int, float, bool, type(None))):
        return obj
    if isinstance(obj, (set, frozenset, list, tuple)):
        return str(sorted([str(x) for x in obj]))
    return str(obj)


def normalize_celltype(name):
    if pd.isna(name) or not isinstance(name, str):
        return None
    norm = name.lower().strip()
    norm = norm.replace('-', ' ')
    norm = re.sub(r'\bstem cells\b', 'stem cell', norm)
    norm = re.sub(r'\bcells\b', 'cell', norm)
    norm = ' '.join(norm.split())
    return norm
