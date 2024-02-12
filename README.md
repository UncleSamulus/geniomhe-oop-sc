# oopsc
Reimplementation of core scanpy plot with Plotly powered visualization.

---

A GENIOMHE Master 1 Python project

## Quick setup

### Installation

```bash
git clone git@github.com:UncleSamulus/geniomhe-oop-sc.git oopsc
cd oopsc
pip install .
```

### Usage

```python
import oopsc.scanpy
```

The above import register a `pli` object in scanpy, as an alternative to a limited subset of scanpy plotting utilities.


To use the implemented dynamical plotting function, you could use:
```python
import scanpy as sc
sc.pli.violin(adata, [...]) 
```

## Milestones

- [ ] Implement scanpy violin plot


## Team members

- Océane Saïbou
- Samuel Ortion 
