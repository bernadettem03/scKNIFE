# Installation

Clone the repository to get the example scripts and graphs. Use Python 3.9.21 and install the tutorial dependencies:

```bash
git clone https://github.com/bernadettem03/scKNIFE.git
cd scKNIFE
conda create -n scknife python=3.9.21 pip
conda activate scknife
python -m pip install -e '.[tutorial]'
```

Check the installation:

```python
import scknife
from importlib.metadata import version
print(version("scknife"))
```

Continue to [Run scKNIFE](tutorials/36721-setup.md#start-with-raw-counts) to prepare the dataset and fit the model.
