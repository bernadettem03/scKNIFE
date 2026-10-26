"""User-facing activity-based cell-type assignment workflow."""
from pathlib import Path
import pandas as pd
from .assignment import compute_marker_anova_f,compute_cluster_marker_means,assign_cell_types

def run_celltype_assignment(adata,activity,marker_nodes,cluster_key="leiden_clusters",
                            methods=("highest_as","zscore"),output_dir=None):
    """Assign cluster labels from marker-node activity without reference labels."""
    if isinstance(activity,(str,Path)): activity=pd.read_parquet(activity)
    if activity.shape[1]!=adata.n_obs: activity=activity.T
    common=adata.obs_names.intersection(activity.columns); adata=adata[common].copy(); activity=activity[common]
    if isinstance(marker_nodes,(str,Path)):
        frame=pd.read_csv(marker_nodes).dropna(subset=["cell_type"]); marker_nodes=sorted(set(frame.cell_type.astype(str)))
    found=[m for m in marker_nodes if m in activity.index]
    if not found: raise ValueError("No marker nodes occur in the activity matrix")
    labels=adata.obs[cluster_key]; clusters=sorted(labels.unique())
    anova=compute_marker_anova_f(found,activity,labels); mappings={}; summaries=[]
    for method in methods:
        means=compute_cluster_marker_means(activity,found,labels,method=method)
        assigned,scores=assign_cell_types(means); col=f"bgnmf_cell_type_{method}"
        adata.obs[col]=labels.map(lambda c:assigned.get(c,"unknown"))
        mapping=pd.DataFrame([{"Cluster":int(c),"Assigned_CellType":assigned.get(c,"unknown"),
                               "Mean_Score":float(scores.get(c,0)),"N_Cells":int((labels==c).sum())} for c in clusters])
        mappings[method]=mapping; summaries.append({"method":method,"n_clusters":len(clusters),"n_assigned_cell_types":adata.obs[col].nunique()})
    summary=pd.DataFrame(summaries)
    if output_dir is not None:
        out=Path(output_dir); out.mkdir(parents=True,exist_ok=True); anova.to_csv(out/'marker_anova_f_stats.csv',index=False); summary.to_csv(out/'assignment_summary.csv',index=False)
        for method,mapping in mappings.items():
            mapping.to_csv(out/f'cluster_celltype_mapping_{method}.csv',index=False)
            adata.obs[[f"bgnmf_cell_type_{method}"]].to_csv(out/f'celltype_assignment_{method}.csv')
        adata.write_h5ad(out/'celltype_assigned.h5ad')
    return adata,mappings,summary
