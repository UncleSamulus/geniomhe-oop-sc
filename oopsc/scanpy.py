"""Add pli to scanpy"""

from . import pli

def patch():
    try:
        import scanpy as sc
    except:
        raise ImportError("Could not patch dynamic plotting API into scanpy")

    setattr(sc, 'pli', pli)

patch()