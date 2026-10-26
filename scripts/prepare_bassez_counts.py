"""Convert Matrix Market counts exported by prepare_bassez_counts.R to AnnData."""
import sys
from pathlib import Path
from scipy.io import mmread
import anndata as ad
import pandas as pd
source=Path(sys.argv[1]);target=Path(sys.argv[2])
x=mmread(source/'counts.mtx').T.tocsr()
genes=(source/'genes.txt').read_text().splitlines()
cells=(source/'cells.txt').read_text().splitlines()
assert x.shape==(len(cells),len(genes))
assert len(set(cells))==len(cells) and len(set(genes))==len(genes)
adata=ad.AnnData(x,obs=pd.DataFrame(index=cells),var=pd.DataFrame(index=genes))
adata.write_h5ad(target)
print(adata.shape)
