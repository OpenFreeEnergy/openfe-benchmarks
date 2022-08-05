import numpy as np
import matplotlib.pylab as plt

from arsenic import plotting, wrangle

fe = wrangle.FEMap('arsenic.csv')

plotting.plot_DDGs(fe.graph, target_name='TYK2', title='Relative Binding Free Energies (TYK2)', filename='./tyk2_relative.png', map_positive=True, figsize=10, color='blue')

experimental_mean_dg = np.asarray([node[1]["exp_DG"] for node in fe.graph.nodes(data=True)]).mean()

plotting.plot_DGs(fe.graph, target_name='TYK2', title='Absolute Binding Free Energies (TYK2)', filename='./tyk2_absolute.png', figsize=10, color='blue')
