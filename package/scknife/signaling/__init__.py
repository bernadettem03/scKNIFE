"""Activity-based ligand-receptor signaling APIs."""
import sys
from pathlib import Path
from . import activity_lr, common

def run_activity_lr(activity_path,h5ad_path,interaction_db,output_dir,permutations=1000,
                    chunk_cells=512,checkpoint_every=50,skip_activity_export=True,
                    group_key="bgnmf_cell_type_highest_as",patient_key="patient_id",
                    condition_key="timepoint",response_key="expansion"):
    """Compute signaling tables with resumable, label-signature-safe permutations."""
    output_dir=Path(output_dir)
    values={"ACTIVITY_PATH":Path(activity_path),"H5AD_PATH":Path(h5ad_path),"CELLCHAT_DB_PATH":Path(interaction_db),
            "DATA_DIR":output_dir,"FIG_DIR":output_dir/"figures","LOG_DIR":output_dir/"logs",
            "GROUP_COLUMN":group_key,"PATIENT_COLUMN":patient_key,
            "SITE_COLUMN":condition_key,"METASTATIC_COLUMN":response_key}
    for module in (common,activity_lr):
        for key,value in values.items(): setattr(module,key,value)
    common.ensure_directories()
    old=sys.argv[:]
    try:
        sys.argv=["activity_lr","--permutations",str(permutations),"--chunk-cells",str(chunk_cells),"--checkpoint-every",str(checkpoint_every)]
        if skip_activity_export: sys.argv.append("--skip-activity-export")
        activity_lr.main()
    finally: sys.argv=old
    return output_dir/"run_manifest.json"

__all__=["run_activity_lr"]
