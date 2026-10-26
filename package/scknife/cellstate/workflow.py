"""User-facing activity-based cell-state assignment workflow."""
from pathlib import Path
import pandas as pd
from .assignment import compute_cluster_state_means,build_state_mapping_table

def run_cellstate_assignment(adata,activity,state_nodes,cluster_key="leiden_clusters",output_dir=None):
    """Assign each cluster to its strongest available state-activity node."""
    if isinstance(activity,(str,Path)): activity=pd.read_parquet(activity)
    if activity.shape[1]!=adata.n_obs: activity=activity.T
    common=adata.obs_names.intersection(activity.columns); adata=adata[common].copy(); activity=activity[common]
    if isinstance(state_nodes,(str,Path)):
        frame=pd.read_csv(state_nodes); state_nodes=sorted(set(frame.cell_state.astype(str)))
    found=[s for s in state_nodes if s in activity.index]
    if not found: raise ValueError("No state nodes occur in the activity matrix")
    for state in found: adata.obs[state]=activity.loc[state]
    assignments,scores=compute_cluster_state_means(adata,found,cluster_col=cluster_key)
    adata.obs["bgnmf_cell_state"]=adata.obs[cluster_key].map(lambda c:assignments[c]["cell_state"])
    mapping=build_state_mapping_table(assignments)
    if output_dir is not None:
        out=Path(output_dir); out.mkdir(parents=True,exist_ok=True); mapping.to_csv(out/'cluster_cellstate_mapping.csv',index=False); scores.to_csv(out/'cluster_cellstate_scores.csv')
        adata.obs[[cluster_key,"bgnmf_cell_state"]].to_csv(out/'cellstate_assignment.csv'); adata.write_h5ad(out/'cellstate_assigned.h5ad')
    return adata,mapping,scores
