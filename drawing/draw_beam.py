# -*- coding: utf-8 -*-
__author__ = 'jakabgabor'

from drawing import _plotting_available, plt


def draw_structure(structure, show=True, deformed=True):

    if _plotting_available:
        for beam in structure.beams:
            # line width is in pixel
            # node size is in pixel
            plt.plot([p.x for p in beam.nodes], [p.y for p in beam.nodes], 'b-', linewidth=2)  # beams
            plt.scatter([p.x for p in beam.nodes], [p.y for p in beam.nodes], marker='s', color='r')  # nodes

        # plot supports
        # 1/scale will be the length of the bars representing the restricted translational degrees of freedom
        scale = 20
        figsize = [plt.rcParams["figure.dpi"] * x for x in plt.rcParams["figure.figsize"]][0]  # width of the fig in pixels
        for k, v in structure.supports.items():
            _aktnode = [x for x in structure.nodes if x.ID == k][0]
            for dof in v:
                if dof == 'ux':  # horizonatal
                    plt.plot([_aktnode.x, _aktnode.x + figsize / scale], [_aktnode.y, _aktnode.y], 'b-', linewidth=4, alpha=0.5)  # a horizontal line
                if dof == 'uy':  # vertical
                    plt.plot([_aktnode.x, _aktnode.x], [_aktnode.y, _aktnode.y - figsize / scale], 'b-', linewidth=4, alpha=0.5)  # a horizontal line
                if dof == 'rotz':  # rotation about Z
                    plt.plot([_aktnode.x, _aktnode.x], [_aktnode.y, _aktnode.y], 'k.', markersize=16, alpha=0.5)  # a point




        # plot loads

        if deformed:
            pass

            # mp = self.midpoint
            # ax = plt.gca()
            # ax.arrow(mp[0], mp[1], self.unitvector[0] * (self.norm / 3.), self.unitvector[1] * (self.norm / 3.),
            #          head_width=self.norm / 30., head_length=self.norm / 15., fc='k', ec='k')
        if show:
            plt.axis('tight')
            plt.axis('equal')
            plt.show()
