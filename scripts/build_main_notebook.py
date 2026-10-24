"""Build the main notebook from the reviewed tutorial pages; --check detects drift."""
import argparse, json, re
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
TARGET=ROOT/'tutorial/scKNIFE_Bassez_Cohort_B_complete_tutorial.ipynb'
PAGES=['setup','annotation','metabolic','go','factors','signaling','validation']
cells=[]
def add(kind,text,**meta):
 text=text.strip()
 if not text:return
 c={'cell_type':kind,'metadata':meta,'source':text.splitlines(keepends=True)}
 if kind=='code':c.update(execution_count=None,outputs=[])
 cells.append(c)
def clean(text):
 text=re.sub(r'^\[(?:Tutorial overview|Previous:).*$', '', text,flags=re.M)
 text=re.sub(r'\{doc\}`([^`]+)`',lambda m:'[Documentation](../docs/'+m[1].replace('../','')+'.md)',text)
 text=re.sub(r'\{download\}`([^<]+)<([^>]+)>`',lambda m:'['+m[1].strip()+'](../docs/'+m[2].replace('../','')+')',text)
 text=re.sub(r'\]\((10000-[^)]+\.md)(#[^)]*)?\)',r'](../docs/tutorials/\1\2)',text)
 text=text.replace('](../downloads/','](../docs/downloads/').replace('](../tutorial-','](../docs/tutorial-')
 return text
add('markdown','''# scKNIFE: complete 10,000-cell Bassez tutorial

Fit scKNIFE on 10,000 Bassez cells, explore clustering, then choose an analysis.

**Get started:** prepare inputs → tune parameters → fit scKNIFE → explore clustering and cell-type retention.

**Analysis tasks:** metabolic tasks, GO analysis, factor interpretation and cell–cell signaling. Each task uses the saved fit and activities; they do not depend on one another.

Use the Python 3.9 environment from the [installation guide](../docs/installation.md) and install the tutorial dependencies:

```bash
python -m pip install -e '.[tutorial]'
```

Follow [input preparation](../docs/tutorial-inputs.md) to obtain and convert the cohort-2 counts, matching metadata and pinned CellChat interaction export. Set `SCKNIFE_COUNTS`, `SCKNIFE_METADATA` and `SCKNIFE_INTERACTION_DB` before launching Jupyter, or set them in the first cell. Input data are not bundled with this notebook.

Open this notebook in Jupyter, or use `python scripts/run_tutorial.py` to save an executed copy. Choose a new output directory for a fresh analysis.''')
add('code','''import os
from pathlib import Path
ROOT = next(p for p in [Path.cwd(), *Path.cwd().parents]
            if (p / "graphs/version_1/unified_graph.graphml").exists())
os.chdir(ROOT)
# Uncomment and set your own input paths if not already set in the environment.
# os.environ["SCKNIFE_COUNTS"] = str(ROOT / "tutorial/inputs/bassez_cohort2_counts.h5ad")
# os.environ["SCKNIFE_METADATA"] = "/absolute/path/1871-BIOKEY_metaData_cohort2_web.csv"
# os.environ["SCKNIFE_INTERACTION_DB"] = str(ROOT / "tutorial/inputs/interactions.csv")
for name in ("SCKNIFE_COUNTS", "SCKNIFE_METADATA", "SCKNIFE_INTERACTION_DB"):
    if not os.environ.get(name) or not Path(os.environ[name]).is_file():
        raise FileNotFoundError(f"Set {name} to the prepared input file; see docs/tutorial-inputs.md")
os.environ.setdefault("SCKNIFE_TUTORIAL_OUTPUT", str(ROOT / "tutorial/outputs_10000"))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(name, "2")
os.environ["MPLBACKEND"] = "Agg"
from IPython.display import Image, display

def show_generated(relative_path):
    path = Path(os.environ["SCKNIFE_TUTORIAL_OUTPUT"]) / relative_path
    if not path.is_file():
        raise FileNotFoundError(f"Run the preceding analysis cell to generate {path.name}")
    display(Image(filename=str(path), width=850))''',tutorial_role='configuration')
for page in PAGES:
 if page=='metabolic':
  add('markdown', '## Analysis tasks\n\nChoose any of the following tasks using the saved fit and activities. Each section loads its own inputs.')
 text=(ROOT/f'docs/tutorials/10000-{page}.md').read_text()
 if page=='setup':
  a=text.index('### Configure the tutorial');b=text.index('## Sample and fit')
  text=text[:a]+text[b:]
 if page=='validation':text=text.split('## Clean-run verification')[0]
 for chunk in re.split(r'(```[^\n]*\n.*?```)',text,flags=re.S):
  if chunk.startswith('```python\n'):
   add('code',chunk[len('```python\n'):-3],tutorial_stage=page)
  elif chunk.startswith('```{image}'):
   filename=chunk.splitlines()[0].split('tutorial-10000/')[1]
   if filename=='organization_umap.png':filename='cellular_organization/'+filename
   if filename.startswith('signaling-seven/'):
    name=Path(filename).name
    mapping={'permutation_fdr.png':'inference/permutation_fdr.png','network_centrality.png':'inference/network_centrality.png','on_minus_pre.png':'inference/on_minus_pre.png','01_sender_receiver.png':'general/01_sender_receiver.png','06_hierarchical_network.png':'general/06_hierarchical_network.png','VEGF_heatmap.png':'pathways/VEGF/03_communication_heatmap.png','VEGF_top_lr.png':'pathways/VEGF/05_top_lr_pairs.png'}
    filename='signaling/'+mapping.get(name,name)
   add('code',f'show_generated({filename!r})',tutorial_role='figure',figure_path=filename)
  elif chunk.startswith('```bash'):continue
  else:add('markdown',clean(chunk))
add('markdown','''## Further resources

See [biological resources and citations](../docs/graph-resources.md) for the resources used in these examples.''')
nb={'cells':cells,'metadata':{'kernelspec':{'display_name':'Python 3 (scKNIFE)','language':'python','name':'python3'},'language_info':{'name':'python','version':'3.9.21'}},'nbformat':4,'nbformat_minor':4}
if '--check' in __import__('sys').argv:
 old=json.loads(TARGET.read_text())
 assert [(c['cell_type'],''.join(c['source'])) for c in old['cells']]==[(c['cell_type'],''.join(c['source'])) for c in cells], 'Notebook differs from reviewed pages; regenerate with scripts/build_main_notebook.py'
 print('Main notebook matches all seven reviewed tutorial pages.')
else:
 TARGET.write_text(json.dumps(nb,indent=1)+'\n')
 print(f'Built {len(cells)} cells, {sum(c["metadata"].get("tutorial_role")=="figure" for c in cells)} generated figures')
