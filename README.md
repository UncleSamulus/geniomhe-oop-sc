# Project Single-Cell

A reimplementation of basic scanpy plotting utility with HoloViews powered dynamical plots.

## Quick setup

### Installation

```bash
git clone git@github.com:UncleSamulus/geniomhe-oop-sc.git
cd geniomhe-oop-sc
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
[...]
sc.pli.violin(adata, [...]) 
```

## Milestones

- [ ] Implement scanpy violin plot


## Team members

- Océane Saïbou
- Samuel Ortion 
