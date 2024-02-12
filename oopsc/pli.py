import scanpy
import plotly

class Plotter:

    def __init__(self):
        pass
    
    def violin(self, adata, keys):
        plots = []
        for key in keys:
            violin = adata.obs.hvplot.violin(y=keys, shared_axes=False)
            plots.append(violin)
        layout = hv.Layout(plots)
        layout.opts(shared_axes=False)
        return layout
    
    def tsne():
        """
        """
    
