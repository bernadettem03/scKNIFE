import os
import json
import time
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from sklearn.decomposition import NMF as SklearnNMF
from sklearn.metrics import calinski_harabasz_score, davies_bouldin_score
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from scknife.tuning.advanced import leiden_silhouette, evaluate_h_graph_clustering
from scknife.utils.helpers import sanity_summary


def plot_loss_curve(loss_history, title, save_path):
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(loss_history, linewidth=1.0)
    ax.set_xlabel('Epoch')
    ax.set_ylabel('Normalized Loss')
    ax.set_title(title, fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, len(loss_history) - 1)
    plt.tight_layout()
    plt.savefig(save_path, dpi=100, bbox_inches='tight')
    plt.close()


def load_best_params(json_path):
    with open(json_path, 'r') as f:
        params = json.load(f)
    print(f"  Loaded best params from {json_path}:")
    for k, v in params.items():
        print(f"    {k} = {v}")
    return params


def save_best_params(params, patient_id, param_name, output_dir):
    out_path = os.path.join(output_dir, f"best_params_{param_name}_{patient_id}.json")
    with open(out_path, 'w') as f:
        json.dump(params, f, indent=2)
    print(f"  Saved best params to {out_path}")
    return out_path


def run_tuning(adata, G, patient_id, params, param_name, param_values,
               seeds=1, n_iter=500, skip_stability=False, output_dir=None,
               resolutions=None, patience=10, defer_stability=False):
    from scknife.graph.laplacian import build_graph_laplacian
    from scknife.graph.mappings import build_Rg_mapping
    from scknife.nmf.initialization import run_nmf, initialize_W
    from scknife.nmf.propagation import propagate_embeddings
    from scknife.nmf.optimization import run_graph_regularized_nmf

    if output_dir is None:
        output_dir = "."
    plots_dir = os.path.join(output_dir, "plots")
    os.makedirs(plots_dir, exist_ok=True)

    K = adata.X.T
    Rg, _, _ = build_Rg_mapping(adata, G)
    A, D, L = build_graph_laplacian(G)

    seed_list = list(range(seeds))
    sweep_values = param_values[param_name]
    total_runs = len(sweep_values) * len(seed_list)
    run_count = 0
    tuning_rows = []
    all_seed_results = {}
    is_k_phase = (param_name == "k")
    prev_sil = None
    prev_stab = None
    prev_recon = None
    consecutive_drops = 0
    early_stopped = False

    if defer_stability:
        skip_stability = True

    print(f"\n  Sweeping {param_name} over {len(sweep_values)} values x {len(seed_list)} seeds = {total_runs} runs\n")

    for val in sweep_values:
        if early_stopped:
            print(f"\n  >>> Early stop triggered after {consecutive_drops} consecutive drops\n")
            break
        p = dict(params)
        p[param_name] = int(val) if param_name == "k" else val
        k = p["k"]
        lam = p["lambda_reg"]
        bw = p["beta_W"]
        bh = p["beta_H"]

        val_str = f"{param_name}={val}"
        seed_results = []
        all_seed_results[val] = seed_results

        # Standard NMF baseline for this k (deterministic, run once)
        nmf_std = SklearnNMF(n_components=k, init='nndsvda', max_iter=1000, random_state=42)
        W_std = nmf_std.fit_transform(K)
        H_std = nmf_std.components_
        X_hat_std = W_std @ H_std
        K_flat = K.flatten().astype(np.float64)
        recon_corr_std_nmf, _ = pearsonr(K_flat, X_hat_std.flatten().astype(np.float64))

        for seed in seed_list:
            run_count += 1
            np.random.seed(seed)
            t_seed_start = time.time()
            print(f"    [{val_str}] seed={seed}  ({run_count}/{total_runs})")

            W_init, H_init, overlap_idx = run_nmf(K, Rg, k=k)
            W_labeled, labeled_idx = initialize_W(W_init, Rg, overlap_idx, Rg.shape[1], k)
            W_prop, _ = propagate_embeddings(W_labeled, L, labeled_idx, Rg.shape[1], k)

            adata_tune = adata.copy()
            W_final, H_final, loss_hist = run_graph_regularized_nmf(
                K, adata_tune, W_prop, Rg, A, D, L,
                k=k, lambda_reg=lam, beta_W=bw, beta_H=bh, n_iter=n_iter
            )

            recon = np.linalg.norm(K - Rg @ W_final @ H_final, 'fro') ** 2
            relative_err = recon / (np.linalg.norm(K, 'fro') ** 2)
            X_hat_scknife = Rg @ W_final @ H_final
            recon_corr_scknife, _ = pearsonr(K_flat, X_hat_scknife.flatten().astype(np.float64))
            sparsity_W = float(np.mean(np.abs(W_final) < 1e-8))
            sparsity_H = float(np.mean(np.abs(H_final) < 1e-8))

            H_cols = [c for c in adata_tune.obs.columns if c.startswith("H_graph_")]
            H_mat = adata_tune.obs[H_cols].values

            sil, sil_res, sil_clusters, labels = leiden_silhouette(H_mat, resolutions=resolutions)
            sil = sil if labels is not None else np.nan

            n_unique = len(np.unique(labels)) if labels is not None else 0
            ch  = calinski_harabasz_score(H_mat, labels) if labels is not None and n_unique > 1 else np.nan
            db  = davies_bouldin_score(H_mat, labels) if labels is not None and n_unique > 1 else np.nan

            sil_str = f"{sil:.4f}" if not np.isnan(sil) else "N/A"
            t_seed = time.time() - t_seed_start
            n_epochs = len(loss_hist)
            print(f"      scKNIFE_corr={recon_corr_scknife:.4f}  StdNMF_corr={recon_corr_std_nmf:.4f}  "
                  f"Sil={sil_str} (res={sil_res})  CH={ch:.1f}  DB={db:.4f}  "
                  f"Loss={loss_hist[-1]:.4f}  Epochs={n_epochs}  [{t_seed:.1f}s]  "
                  f"lambda={lam}  beta_W={bw}  beta_H={bh}")

            plot_title = f"{param_name}={val} seed={seed} (lambda={lam}, beta_W={bw}, beta_H={bh})"
            plot_fname = f"loss_{param_name}_{val}_seed{seed}.png"
            plot_loss_curve(loss_hist, plot_title, os.path.join(plots_dir, plot_fname))

            row = {
                "seed": seed, "recon_corr_scknife": recon_corr_scknife,
                "recon_corr_std_nmf": recon_corr_std_nmf,
                "silhouette": sil,
                "calinski_harabasz": ch, "davies_bouldin": db,
                "recon_error": recon, "relative_error": relative_err,
                "sparsity_W": sparsity_W, "sparsity_H": sparsity_H,
                "loss_final": loss_hist[-1], "n_epochs": n_epochs, "time_s": t_seed,
                "lambda_reg": lam, "beta_W": bw, "beta_H": bh,
                "W": W_final, "H": H_final,
                param_name: int(val) if param_name == "k" else val,
            }
            seed_results.append(row)

        if not skip_stability:
            Ws = [r["W"] for r in seed_results]
            corr_list = []
            for i in range(len(seed_list)):
                for j in range(i + 1, len(seed_list)):
                    top_k = min(10, k)
                    c_vals = []
                    for c in range(top_k):
                        try:
                            c_vals.append(np.corrcoef(Ws[i][:, c], Ws[j][:, c])[0, 1])
                        except Exception:
                            c_vals.append(0)
                    corr_list.append(np.mean(c_vals))
            stability = float(np.mean(corr_list)) if corr_list else 0.0
        else:
            stability = np.nan

        metrics = ["recon_corr_scknife", "recon_corr_std_nmf", "silhouette",
                   "calinski_harabasz", "davies_bouldin",
                   "recon_error", "relative_error", "sparsity_W", "sparsity_H",
                   "loss_final", "n_epochs", "time_s"]
        avg = {m: float(np.nanmean([r[m] for r in seed_results])) for m in metrics}
        avg[param_name] = int(val) if param_name == "k" else val
        avg["stability"] = stability
        avg["lambda_reg"] = lam
        avg["beta_W"] = bw
        avg["beta_H"] = bh
        tuning_rows.append(avg)

        metric_line = (f"CORR_SCKNIFE={avg['recon_corr_scknife']:.3f}  "
                       f"CORR_STD_NMF={avg['recon_corr_std_nmf']:.3f}  "
                       f"SIL={avg['silhouette']:.3f}  "
                       f"Loss={avg['loss_final']:.4f}  "
                       f"Epochs={int(avg['n_epochs'])}  "
                       f"Time={avg['time_s']:.1f}s")
        stab_str = f"  STAB={stability:.3f}" if not defer_stability and not np.isnan(stability) else ""
        print(f"    -> [{val_str}] {metric_line}{stab_str}\n")

        # Early stopping for k-phase only
        if is_k_phase and not np.isnan(avg["silhouette"]):
            if defer_stability:
                if prev_sil is not None and prev_recon is not None:
                    recon_dropped = avg["recon_corr_scknife"] < prev_recon
                    sil_dropped = avg["silhouette"] < prev_sil
                    if recon_dropped and sil_dropped:
                        consecutive_drops += 1
                        print(f"      Early stop check: drop {consecutive_drops}/{patience} (recon={recon_dropped} sil={sil_dropped})")
                        if consecutive_drops >= patience:
                            tuning_rows = tuning_rows[:-consecutive_drops]
                            for v in list(all_seed_results.keys())[-consecutive_drops:]:
                                del all_seed_results[v]
                            early_stopped = True
                    else:
                        consecutive_drops = 0
            else:
                if not np.isnan(stability) and prev_sil is not None and prev_stab is not None and prev_recon is not None:
                    recon_dropped = avg["recon_corr_scknife"] < prev_recon
                    sil_dropped = avg["silhouette"] < prev_sil
                    stab_dropped = stability < prev_stab
                    if recon_dropped and sil_dropped and stab_dropped:
                        consecutive_drops += 1
                        print(f"      Early stop check: drop {consecutive_drops}/{patience} (recon={recon_dropped} sil={sil_dropped} stab={stab_dropped})")
                        if consecutive_drops >= patience:
                            tuning_rows = tuning_rows[:-consecutive_drops]
                            for v in list(all_seed_results.keys())[-consecutive_drops:]:
                                del all_seed_results[v]
                            early_stopped = True
                    else:
                        consecutive_drops = 0
            prev_sil = avg["silhouette"]
            prev_stab = stability if not defer_stability else np.nan
            prev_recon = avg["recon_corr_scknife"]

        for r in seed_results:
            r.pop("W", None)
            r.pop("H", None)

    all_raw = []
    for val in sweep_values:
        for sr in all_seed_results[val]:
            all_raw.append(sr)
    raw_df = pd.DataFrame(all_raw)
    raw_csv = os.path.join(output_dir, f"tuning_raw_{param_name}_{patient_id}.csv")
    raw_df.to_csv(raw_csv, index=False)

    summary_df = pd.DataFrame(tuning_rows).sort_values(by=param_name)
    summary_csv = os.path.join(output_dir, f"tuning_summary_{param_name}_{patient_id}.csv")
    summary_df.to_csv(summary_csv, index=False)

    if param_name == "beta_W":
        lower, upper = 0.35, 0.55
        candidates = summary_df[
            (summary_df["sparsity_W"] >= lower) & (summary_df["sparsity_W"] <= upper)
        ].copy()
        if candidates.empty:
            print(f"  WARNING: No beta_W value gives sparsity_W in [{lower}, {upper}]. "
                  f"Picking closest to range.")
            df_copy = summary_df.copy()
            df_copy["sparsity_dist"] = df_copy["sparsity_W"].apply(
                lambda x: min(abs(x - lower), abs(x - upper))
            )
            best_row = df_copy.loc[df_copy["sparsity_dist"].idxmin()]
        else:
            best_row = candidates.loc[candidates["silhouette"].idxmax()]
            print(f"  beta_W sparsity filter [{lower}, {upper}]: "
                  f"{len(candidates)}/{len(summary_df)} values pass. "
                  f"Selected beta_W={best_row['beta_W']:.4f} "
                  f"(sparsity_W={best_row['sparsity_W']:.3f}, silhouette={best_row['silhouette']:.4f})")
    else:
        best_row = summary_df.loc[summary_df["silhouette"].idxmax()]
    best_val = int(best_row[param_name]) if param_name == "k" else best_row[param_name]
    print(f"  \u2713 Raw results: {raw_csv}")
    print(f"  \u2713 Summary:     {summary_csv}")
    selection_method = "sparsity range [0.35, 0.55] + silhouette" if param_name == "beta_W" else "silhouette"
    print(f"  \u2192 Optimal {param_name} = {best_val} (by {selection_method})")
    print(f"     Best scKNIFE_corr = {best_row.get('recon_corr_scknife', float('nan')):.4f}")
    print(f"     Best StdNMF_corr  = {best_row.get('recon_corr_std_nmf', float('nan')):.4f}")
    print(f"     Best sparsity_W   = {best_row.get('sparsity_W', float('nan')):.4f}")
    print(f"     Best sparsity_H   = {best_row.get('sparsity_H', float('nan')):.4f}")

    # ── Summary table ──
    print(f"\n")
    print(f"  " + "=" * 100)
    print(f"  Tuning Summary: {param_name}")
    print(f"  " + "=" * 100)
    headers = [param_name, "scKNIFE_corr", "StdNMF_corr", "Silhouette", "sparsity_W", "sparsity_H", "Loss", "Epochs", "Time(s)", "lambda_reg", "beta_W", "beta_H"]
    col_w = max(len(h) for h in headers) + 2
    header_line = "  " + " | ".join(h.rjust(col_w) for h in headers)
    print(header_line)
    print("  " + "-" * len(header_line))
    for r in tuning_rows:
        param_val = int(r[param_name]) if param_name == "k" else r[param_name]
        is_best = (param_val == best_val)
        marker = ">>" if is_best else "  "
        vals = [
            str(param_val),
            f"{r['recon_corr_scknife']:.4f}",
            f"{r['recon_corr_std_nmf']:.4f}",
            f"{r['silhouette']:.4f}",
            f"{r['sparsity_W']:.4f}",
            f"{r['sparsity_H']:.4f}",
            f"{r['loss_final']:.4f}",
            str(int(r['n_epochs'])),
            f"{r['time_s']:.1f}",
            f"{r['lambda_reg']}",
            f"{r['beta_W']}",
            f"{r['beta_H']}",
        ]
        print(marker + " | ".join(v.rjust(col_w) for v in vals))
    print(f"  " + "=" * 100)

    sanity_summary("Tuning", [
        ("param", True, f"{param_name}"),
        ("values", True, f"{len(sweep_values)} x {len(seed_list)} seeds = {total_runs} runs"),
        ("optimal", True, f"{param_name}={best_val}"),
    ])

    # Post-sweep stability: re-run best value with 3 seeds
    if defer_stability and not np.isnan(best_row.get("silhouette", np.nan)):
        stab_seeds = 3
        p_best = dict(params)
        p_best[param_name] = int(best_val) if param_name == "k" else best_val
        k_best = p_best["k"]
        lam_best = p_best["lambda_reg"]
        bw_best = p_best["beta_W"]
        bh_best = p_best["beta_H"]

        print(f"\n  Computing stability for best {param_name}={best_val} "
              f"(re-running {stab_seeds} seeds)...")
        Ws = []
        for s in range(stab_seeds):
            np.random.seed(s)
            W_i, H_i, idx_i = run_nmf(K, Rg, k=k_best)
            W_l, lab_i = initialize_W(W_i, Rg, idx_i, Rg.shape[1], k_best)
            W_p, _ = propagate_embeddings(W_l, L, lab_i, Rg.shape[1], k_best)
            adata_tmp = adata.copy()
            W_f, H_f, _ = run_graph_regularized_nmf(
                K, adata_tmp, W_p, Rg, A, D, L,
                k=k_best, lambda_reg=lam_best,
                beta_W=bw_best, beta_H=bh_best, n_iter=n_iter)
            Ws.append(W_f)

        corr_list = []
        for i in range(len(Ws)):
            for j in range(i + 1, len(Ws)):
                top_k = min(10, k_best)
                c_vals = []
                for c in range(top_k):
                    try:
                        c_vals.append(np.corrcoef(Ws[i][:, c], Ws[j][:, c])[0, 1])
                    except Exception:
                        c_vals.append(0)
                corr_list.append(np.mean(c_vals))
        stability = float(np.mean(corr_list)) if corr_list else 0.0
        print(f"  Stability for {param_name}={best_val}: {stability:.4f}")

        # Update summary CSV with stability
        summary_df["stability"] = stability
        summary_df.to_csv(summary_csv, index=False)

    return best_val, raw_df, summary_df


def run_joint_k_lambda_grid(adata, G, candidate_k, lambda_values, default_params,
                             n_neighbors=20, metric="cosine", resolution=0.5,
                             sample_size=2000):
    """
    Joint grid search over k and lambda_reg values on H_graph directly.

    For each (k, lambda) combination:
      1. Run full scKNIFE pipeline on a copy of adata.
      2. Compute activity scores.
      3. Extract H_graph from obs.
      4. Evaluate clustering quality via evaluate_h_graph_clustering
         (neighbors on X_H_raw + Leiden + silhouette, NO PCA).
      5. Capture reconstruction correlation from pipeline metrics + loss history.
      6. Track best combination by silhouette.

    Parameters
    ----------
    adata : AnnData
        Preprocessed (subsampled) data for tuning.
    G : nx.Graph
        Unified knowledge graph.
    candidate_k : list of int
        k values to evaluate.
    lambda_values : list of float
        lambda_reg values to evaluate.
    default_params : dict
        Must contain "beta_W", "beta_H", "n_iter".
    n_neighbors : int
        Neighbors for H_graph graph construction (default 20).
    metric : str
        Distance metric (default "cosine").
    resolution : float
        Leiden resolution (default 0.5).
    sample_size : int
        Silhouette sample size (default 2000).

    Returns
    -------
    results_df : pd.DataFrame
        Silhouette, n_clusters, recon_corr for each (k, lambda) combination.
    best_k : int
        Best k by silhouette.
    best_lambda : float
        Best lambda_reg by silhouette.
    """
    from scknife import run_full_pipeline, compute_activity_scores

    rows = []
    best_sil = -1.0
    best_k = None
    best_lambda = None
    total = len(candidate_k) * len(lambda_values)
    count = 0

    plots_dir = os.path.join(os.getcwd(), "hyperparameter_tuning", "plots")
    os.makedirs(plots_dir, exist_ok=True)

    print(f"\n  Joint k x lambda grid: {len(candidate_k)} k x {len(lambda_values)} lambda = {total} runs\n"
          f"  Evaluation: neighbors(n={n_neighbors}, metric={metric}) + Leiden(r={resolution}) + silhouette\n")

    for k in candidate_k:
        for lam in lambda_values:
            count += 1
            lam_str = f"{lam}" if isinstance(lam, int) else f"{lam:.1e}"
            t_start = time.time()
            print(f"    [{count}/{total}] k={k}, lambda_reg={lam_str} ...", end=" ", flush=True)

            adata_k = adata.copy()
            adata_k = run_full_pipeline(
                adata=adata_k, G=G, k=k, lambda_reg=lam,
                beta_W=default_params.get("beta_W", 0.01),
                beta_H=default_params.get("beta_H", 0.1),
                n_iter=default_params.get("n_iter", 500),
                verbose=False
            )
            adata_k = compute_activity_scores(adata_k)

            t_elapsed = time.time() - t_start

            # Extract reconstruction correlation
            recon_metrics = adata_k.uns.get("recon_metrics", {})
            recon_corr_scknife = recon_metrics.get("scknife_corr", float("nan"))
            recon_corr_std_nmf = recon_metrics.get("std_nmf_corr", float("nan"))

            # Extract loss history for this run
            loss_hist = adata_k.uns.get("loss_history", [])
            final_loss = loss_hist[-1] if loss_hist else float("nan")
            n_epochs = len(loss_hist)

            # Save loss curve plot
            lam_str_fmt = f"{lam}" if isinstance(lam, int) else f"{lam:.1e}"
            plot_fname = f"loss_joint_k{k}_lambda{lam_str_fmt}.png"
            plot_loss_curve(loss_hist, f"Joint k={k}, lambda_reg={lam_str_fmt} (beta_W={default_params.get('beta_W',0.01)}, beta_H={default_params.get('beta_H',0.1)})", os.path.join(plots_dir, plot_fname))

            # Extract H_graph and evaluate clustering
            h_cols = [c for c in adata_k.obs.columns if c.startswith("H_graph_")]
            if not h_cols:
                print("NO H_GRAPH — skipping")
                continue
            H_mat = adata_k.obs[h_cols].to_numpy()

            cluster_metrics = evaluate_h_graph_clustering(
                H_mat, n_neighbors=n_neighbors, metric=metric,
                resolution=resolution, sample_size=sample_size
            )
            cluster_metrics["k"] = int(k)
            cluster_metrics["lambda_reg"] = float(lam)
            cluster_metrics["recon_corr_scknife"] = recon_corr_scknife
            cluster_metrics["recon_corr_std_nmf"] = recon_corr_std_nmf
            cluster_metrics["loss_final"] = final_loss
            cluster_metrics["n_epochs"] = n_epochs
            cluster_metrics["time_s"] = t_elapsed
            rows.append(cluster_metrics)

            print(f"corr_scknife={recon_corr_scknife:.4f}  corr_std_nmf={recon_corr_std_nmf:.4f}  "
                  f"sil={cluster_metrics['silhouette']:.4f}  n_clust={cluster_metrics['n_clusters']}  "
                  f"loss={final_loss:.4f}  epochs={n_epochs}  [{t_elapsed:.1f}s]")

            if cluster_metrics["silhouette"] > best_sil:
                best_sil = cluster_metrics["silhouette"]
                best_k = int(k)
                best_lambda = float(lam)

    results_df = pd.DataFrame(rows)
    best_row = results_df.loc[results_df["silhouette"].idxmax()] if len(results_df) > 0 else {}
    best_scknife_corr = best_row.get('recon_corr_scknife', float('nan'))
    best_std_nmf_corr = best_row.get('recon_corr_std_nmf', float('nan'))
    best_loss = best_row.get('loss_final', float('nan'))
    best_epochs = best_row.get('n_epochs', float('nan'))
    best_time = best_row.get('time_s', float('nan'))

    # ── Summary table ──
    print(f"\n")
    print(f"  " + "=" * 130)
    print(f"  Joint k x lambda_reg Grid Summary")
    print(f"  " + "=" * 130)
    headers = ["k", "lambda_reg", "scKNIFE_corr", "StdNMF_corr", "Silhouette", "n_clusters", "Loss", "Epochs", "Time(s)"]
    col_w = max(len(h) for h in headers) + 2
    header_line = "  " + " | ".join(h.rjust(col_w) for h in headers)
    print(header_line)
    print("  " + "-" * len(header_line))
    for _, r in results_df.iterrows():
        is_best = (r['k'] == best_k and r['lambda_reg'] == best_lambda)
        marker = ">>" if is_best else "  "
        vals = [
            str(int(r['k'])),
            f"{r['lambda_reg']}",
            f"{r['recon_corr_scknife']:.4f}" if not np.isnan(r['recon_corr_scknife']) else "nan",
            f"{r['recon_corr_std_nmf']:.4f}" if not np.isnan(r['recon_corr_std_nmf']) else "nan",
            f"{r['silhouette']:.4f}",
            str(int(r['n_clusters'])),
            f"{r['loss_final']:.4f}" if not np.isnan(r['loss_final']) else "nan",
            str(int(r['n_epochs'])),
            f"{r['time_s']:.1f}",
        ]
        print(marker + " | ".join(v.rjust(col_w) for v in vals))
    print(f"  " + "=" * 130)

    print(f"\n  >>> Best joint (k, lambda_reg) = ({best_k}, {best_lambda})  "
          f"sil={best_sil:.4f}  "
          f"corr_scknife={best_scknife_corr:.4f}  "
          f"loss={best_loss:.4f}  "
          f"epochs={best_epochs}  "
          f"time={best_time:.1f}s <<<\n")
    return results_df, best_k, best_lambda
