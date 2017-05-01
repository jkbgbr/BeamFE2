# -*- coding: utf-8 -*-
__author__ = 'jakabgabor'

try:
    import matplotlib.pyplot as plt
    _plotting_available = True
    # plt.xkcd()
except ImportError:
    _plotting_available = False
    plt = None
    print('No plotting functions available, matplotlib is missing')
