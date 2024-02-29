# `oopsc`
Reimplementation of core `scanpy` plot with `plotly` powered visualization.

---

A GENIOMHE Master 1 Python project

## Quick setup

### Installation

```bash
git clone git@github.com:UncleSamulus/geniomhe-oop-sc.git oopsc
cd oopsc
python3 -m venv .venv/oopsc
source .venv/oopsc/bin/activate
pip install -e .
```

### Usage

```python
import oopsc.scanpy
```

The above import registers a `pli` object in `scanpy`, as an alternative to `scanpy.pl` for a limited subset of `scanpy` plotting utilities.


To use the implemented dynamical plotting function, you could use:
```python
import scanpy as sc
# Then for example
sc.pli.violin(adata, [...]) 
```

## Team members

- Océane Saïbou
- Samuel Ortion 
