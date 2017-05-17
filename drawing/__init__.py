# -*- coding: utf-8 -*-

try:
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon
    from matplotlib.collections import PatchCollection
    _plotting_available = True
    # plt.xkcd()
except ImportError:
    _plotting_available = False
    plt = None
    print('No plotting functions available, matplotlib is missing')
