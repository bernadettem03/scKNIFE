"""Prevent the documented 10,000-cell commands drifting from executable code."""
from pathlib import Path
import ast
import re
ROOT=Path(__file__).resolve().parents[1]
PAGES={'setup':'','annotation':'_organization','metabolic':'_metabolic','go':'_go','factors':'_factors','signaling':'_signaling','validation':'_validation'}
for page,suffix in PAGES.items():
 text=(ROOT/f'docs/tutorials/10000-{page}.md').read_text()
 blocks=re.findall(r'```python\n(.*?)```',text,re.S)
 if page=='validation':blocks=blocks[:1]
 documented=ast.dump(ast.parse('\n\n'.join(blocks)),include_attributes=False)
 script=ROOT/f'scripts/run_tutorial_10000{suffix}.py'
 assert script.exists(),script
 assert ast.dump(ast.parse(script.read_text()),include_attributes=False)==documented,f'Code drift: {page}'
quick=(ROOT/'docs/quickstart.md').read_text()
for suffix in PAGES.values():
 assert f'python scripts/run_tutorial_10000{suffix}.py' in quick,f'Missing quick-start stage: {suffix}'
setup=(ROOT/'scripts/run_tutorial_10000.py').read_text()
assert 'run_celltype_assignment' not in setup and 'paired_task_comparisons' not in setup
print('All seven tutorial stages match their documented code and quick-start commands.')

import subprocess, sys
subprocess.run([sys.executable, str(ROOT / "scripts/build_main_notebook.py"), "--check"], check=True)

assert (ROOT / "package/scknife/tutorial_helpers.py").read_bytes() == (ROOT / "docs/downloads/factor_selection_helpers.py").read_bytes(), "Factor helper download drift"
assert (ROOT / "scripts/run_tutorial_10000_factors.py").read_bytes() == (ROOT / "docs/downloads/factor_interpretation.py").read_bytes(), "Factor script download drift"
assert (ROOT / "tutorial/scKNIFE_Bassez_Cohort_B_complete_tutorial.ipynb").read_bytes() == (ROOT / "docs/downloads/scKNIFE_10000_tutorial.ipynb").read_bytes(), "Notebook download drift"
assert "run_tutorial_10000_validation.py" not in quick.split("## Optional:")[0], "Projection must stay optional in quick start"

assert (ROOT / "scripts/run_tutorial_10000_metabolic.py").read_bytes() == (ROOT / "docs/downloads/metabolic_walkthrough.py").read_bytes(), "Metabolic script download drift"

assert (ROOT / "scripts/run_tutorial_10000_go.py").read_bytes() == (ROOT / "docs/downloads/go_differential.py").read_bytes(), "GO script download drift"
assert (ROOT / "package/scknife/tutorial_differential.py").read_bytes() == (ROOT / "docs/downloads/differential_activity_helpers.py").read_bytes(), "Differential helper download drift"
