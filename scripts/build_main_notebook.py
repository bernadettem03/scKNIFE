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
add('markdown','''# scKNIFE: 10,000-cell tutorial

Fit scKNIFE on Bassez cohort-2 cells, explore clustering, then choose a biological analysis.

## Getting started

Use the Python 3.9 environment from the [installation guide](../docs/installation.md):

```bash
python -m pip install -e '.[tutorial]'
```

Prepare the [input files](../docs/tutorial-inputs.md) and set `SCKNIFE_COUNTS` and `SCKNIFE_METADATA` in your environment or in the first cell. Choose a new output directory for a fresh fit.''')
add('code','''import os
from pathlib import Path
ROOT = next(p for p in [Path.cwd(), *Path.cwd().parents]
            if (p / "graphs/version_1/unified_graph.graphml").exists())
os.chdir(ROOT)
# Uncomment and set your own input paths if not already set in the environment.
# os.environ["SCKNIFE_COUNTS"] = str(ROOT / "tutorial/inputs/bassez_cohort2_counts.h5ad")
# os.environ["SCKNIFE_METADATA"] = "/absolute/path/1871-BIOKEY_metaData_cohort2_web.csv"
for name in ("SCKNIFE_COUNTS", "SCKNIFE_METADATA"):
    if not os.environ.get(name) or not Path(os.environ[name]).is_file():
        raise FileNotFoundError(f"Set {name} to the prepared input file; see docs/tutorial-inputs.md")
os.environ.setdefault("SCKNIFE_TUTORIAL_OUTPUT", str(ROOT / "tutorial/outputs_10000"))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(name, "2")
os.environ["MPLBACKEND"] = "Agg"
from IPython.display import Image, display

def show_generated(relative_path, width=850):
    path = Path(os.environ["SCKNIFE_TUTORIAL_OUTPUT"]) / relative_path
    if not path.is_file():
        raise FileNotFoundError(f"Run the preceding analysis cell to generate {path.name}")
    display(Image(filename=str(path), width=width))''',tutorial_role='configuration')
for page in PAGES:
 if page=='metabolic':
  add('markdown', '## Analysis tasks\n\nChoose any task below using the saved fit and activities.')
 text=(ROOT/f'docs/tutorials/10000-{page}.md').read_text()
 text=re.sub(r'\A---\n.*?\n---\n+', '', text, flags=re.S)
 # Keep website navigation and repeated prerequisites out of the notebook.
 text=re.sub(r'^\[(?:Tutorial overview|Previous:|Run scKNIFE)\].*$', '', text, flags=re.M)
 if page=='setup':
  text='# Run scKNIFE\n\n'+text[text.index('## 1. Select the 10,000-cell sample'):]
  text=text.split('## Continue with the fitted outputs')[0]
 elif page=='annotation':
  text=text.replace('## Compare the learned organization with source labels\n', '')
  text=re.sub(r' \[Choose an analysis\].*?\n', '\n', text)
 elif page in ('metabolic','go','factors'):
  text=re.sub(r'^\*\*Input:\*\*.*\n', '', text, flags=re.M)
 elif page=='signaling':
  text=re.sub(r'^\*\*Input:\*\*.*\n', '', text, flags=re.M)
 if page=='validation':
  text=text.split('## Clean-run verification')[0]
  text=re.sub(r'^\*\*Input:\*\*.*\n', '', text, flags=re.M)
  text=text.replace('## Rank sensitivity and projection code\n', '')
 # One title, two main sections, then task headings and their code/figure topics.
 titles={'metabolic':'Metabolic tasks','go':'GO analysis'}
 if page in titles:
  text=re.sub(r'^# .*$', '# '+titles[page], text, count=1, flags=re.M)
 renames={
  'metabolic':{'Load activities and choose a comparison':'Load metabolic activities and choose a comparison','Save results':'Save metabolic results'},
  'go':{'Load activities and choose a comparison':'Load GO activities and choose a comparison','2. Examine one cell type':'2. Examine GO activity in one cell type','Save results':'Save GO results'},
  'factors':{'Outputs':'Factor outputs'},
  'signaling':{'Outputs':'Signaling outputs'},
 }
 for old,new in renames.get(page,{}).items():
  text=text.replace('## '+old+'\n','## '+new+'\n')
 offset=1 if page=='validation' else 2
 for chunk in re.split(r'(```[^\n]*\n.*?```)',text,flags=re.S):
  if chunk.startswith('```python\n'):
   add('code',chunk[len('```python\n'):-3],tutorial_stage=page)
  elif chunk.startswith('```{image}'):
   filename=chunk.splitlines()[0].split('tutorial-10000/')[1]
   if filename=='organization_umap.png':filename='cellular_organization/'+filename
   if filename.startswith('signaling-seven/'):
    name=Path(filename).name
    mapping={'permutation_fdr.png':'inference/permutation_fdr.png','network_centrality.png':'inference/network_centrality.png','01_sender_receiver.png':'general/01_sender_receiver.png','06_hierarchical_network.png':'general/06_hierarchical_network.png','VEGF_heatmap.png':'pathways/VEGF/03_communication_heatmap.png','VEGF_top_lr.png':'pathways/VEGF/05_top_lr_pairs.png'}
    filename='signaling/'+mapping.get(name,name)
   width = 380 if filename == 'metabolic_walkthrough/metabolic_task_distribution.png' else 680
   display_code = f'show_generated({filename!r}, width={width})' if page == 'metabolic' else f'show_generated({filename!r})'
   add('code',display_code,tutorial_role='figure',figure_path=filename)
  elif chunk.startswith('```bash'):continue
  else:
   chunk=re.sub(r'^(#{1,4}) ', lambda m: '#'*(len(m[1])+offset)+' ', chunk, flags=re.M)
   add('markdown',re.sub(r'\n{3,}', '\n\n', clean(chunk)))
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
