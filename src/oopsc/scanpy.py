"""Add pli to scanpy"""

from .pli import Plotter

def patch():
    try:
        import scanpy as sc
    except:
        raise ImportError("Could not patch dynamic plotting API into scanpy")

    pli = Plotter()
    setattr(sc, 'pli', pli)

patch()