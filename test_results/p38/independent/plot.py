import numpy as np
import matplotlib.pylab as plt

from arsenic import plotting, wrangle

fe = wrangle.FEMap('arsenic.csv')

plotting.plot_DDGs(fe.graph, target_name='P38', title='Relative Binding Free Energies (P38 independent windows)', filename='./p38_independent_relative.png', map_positive=True, figsize=10, color='blue')

experimental_mean_dg = np.asarray([node[1]["exp_DG"] for node in fe.graph.nodes(data=True)]).mean()

plotting.plot_DGs(fe.graph, target_name='P38', title='Absolute Binding Free Energies (P38 indepenendent windows)', filename='./p38_independent_absolute.png', figsize=10, color='blue')
