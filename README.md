# `oopsc`
Reimplementation of core `scanpy` plot with `plotly` powered visualization.

---

A GENIOMHE Master 1 Python project

## Quick setup

### Installation

```bash
# First clone the repository
git clone https://github.com/samuelortion/geniomhe-oop-sc.git oopsc
cd oopsc
# To install oopsc and its dependencies
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

## Showcase

We develop a dashboard interface with dash.

First, install `dash`
```bash
pip install -e .'showcase'
```

Then, run the dash application:
```bash
python3 showcase/main.py
```

## Team members

- Océane Saïbou
- Samuel Ortion 
