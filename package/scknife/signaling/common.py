#!/usr/bin/env python3
"""Configuration shared by the activity-based ligand-receptor workflow."""
from __future__ import annotations
import hashlib, json, os, platform, sys
from pathlib import Path
import numpy as np

DATA_DIR=Path(os.environ.get("SCKNIFE_SIGNALING_OUTPUT","scknife_signaling_output")).resolve()
FIG_DIR=Path(os.environ.get("SCKNIFE_SIGNALING_FIGURES",DATA_DIR/"figures")).resolve()
LOG_DIR=Path(os.environ.get("SCKNIFE_SIGNALING_LOGS",DATA_DIR/"logs")).resolve()
ACTIVITY_PATH=Path(os.environ.get("SCKNIFE_ACTIVITY_PATH","activity_scores.parquet")).resolve()
H5AD_PATH=Path(os.environ.get("SCKNIFE_H5AD_PATH","scknife_results.h5ad")).resolve()
CELLCHAT_DB_PATH=Path(os.environ.get("SCKNIFE_INTERACTION_DB","interactions.csv")).resolve()
GROUP_COLUMN="bgnmf_cell_type_highest_as"; PATIENT_COLUMN="patient_id"; SITE_COLUMN="timepoint"; METASTATIC_COLUMN="expansion"
RANDOM_SEED=42; N_PERMUTATIONS=1000; FDR_THRESHOLD=0.05
DB_ANNOTATIONS=("Secreted Signaling","ECM-Receptor","Cell-Cell Contact")
PATHWAY_ORDER=("TGFb","CXCL","VEGF","WNT","NOTCH","IFN","TNF")
CELLTYPE_STATE_ORDER=["t memory cell","mammary epithelial cell","dendritic cell","pancreatic stellate cell","endothelial cell"]
H_COLORS={"t memory cell":"#7f7f7f","mammary epithelial cell":"#1f77b4","dendritic cell":"#e377c2","pancreatic stellate cell":"#8c564b","endothelial cell":"#aec7e8"}
def ensure_directories():
 for p in (DATA_DIR,FIG_DIR,LOG_DIR): p.mkdir(parents=True,exist_ok=True)
def write_json(path,payload):
 path.parent.mkdir(parents=True,exist_ok=True); path.write_text(json.dumps(payload,indent=2,sort_keys=True,default=str)+"\n")
def software_versions():
 import anndata,matplotlib,networkx,pandas,pyarrow,scipy,seaborn,statsmodels
 return {"python":sys.version,"platform":platform.platform(),"numpy":np.__version__,"pandas":pandas.__version__,"scipy":scipy.__version__,"statsmodels":statsmodels.__version__,"networkx":networkx.__version__,"matplotlib":matplotlib.__version__,"seaborn":seaborn.__version__,"anndata":anndata.__version__,"pyarrow":pyarrow.__version__,"cnsplots":"unavailable; matplotlib/seaborn fallback"}
def atomic_npz(path,**arrays):
 tmp=path.with_name(path.name+".tmp.npz"); np.savez_compressed(tmp,**arrays); os.replace(tmp,path)
