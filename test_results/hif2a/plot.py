import numpy as np
import matplotlib.pylab as plt

from arsenic import plotting, wrangle

fe = wrangle.FEMap('arsenic.csv')

plotting.plot_DDGs(fe.graph, target_name='HIF2A', title='Relative Binding Free Energies (HIF2A)', filename='./hif2a_relative.png', symmetrise=True, figsize=10, color='blue')

experimental_mean_dg = np.asarray([node[1]["exp_DG"] for node in fe.graph.nodes(data=True)]).mean()

plotting.plot_DGs(fe.graph, target_name='HIF2A', title='Absolute Binding Free Energies (HIF2A)', filename='./hif2a_absolute.png', figsize=10, color='blue')
