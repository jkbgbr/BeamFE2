# -*- coding: utf-8 -*-
__author__ = 'jakabgabor'

from drawing import _plotting_available, plt
import itertools


def draw_structure(structure, show=True, deformed=True):

    if _plotting_available:
        for beam in structure.beams:
            # line width is in pixel
            plt.plot([p.x for p in beam.nodes], [p.y for p in beam.nodes], 'b-', linewidth=3, alpha=0.5, zorder=1)  # beams
            plt.scatter([p.x for p in beam.nodes], [p.y for p in beam.nodes], marker='s', color='r', alpha=0.5, s=50, zorder=2)  # nodes



        # plot supports
        # 1/scale will be the length of the bars representing the restricted translational degrees of freedom
        scale = 20
        figsize = [plt.rcParams["figure.dpi"] * x for x in plt.rcParams["figure.figsize"]]  # width of the fig in pixels
        w_figsize = figsize[0]
        h_figsize = figsize[1]
        for k, v in structure.supports.items():
            _aktnode = [x for x in structure.nodes if x.ID == k][0]
            for dof in v:
                if dof == 'ux':  # horizonatal
                    plt.plot([_aktnode.x, _aktnode.x + w_figsize / scale], [_aktnode.y, _aktnode.y], 'g-', linewidth=4, zorder=3)  # a horizontal line
                if dof == 'uy':  # vertical
                    plt.plot([_aktnode.x, _aktnode.x], [_aktnode.y, _aktnode.y - w_figsize / scale], 'g-', linewidth=4, zorder=3)  # a horizontal line
                if dof == 'rotz':  # rotation about Z
                    plt.plot([_aktnode.x, _aktnode.x], [_aktnode.y, _aktnode.y], 'g.', markersize=16, zorder=3)  # a point

        # plot loads


        # plot the deformed shape
        if deformed:
            fig = plt.gca()
            # length of the longes beam - this will be the base for the scaling
            _long = sorted(structure.beams, key=lambda x: x.l)[-1].l

            # displacements
            dre = structure.resulting_displacement
            # the scaling factor, based on the larges displacement and the length of the longest element
            _scale = (_long / 10.) / max(dre)

            for beam in structure.beams:
                # line width is in pixel
                dxs = beam.displacement_component(component='ux')
                dys = beam.displacement_component(component='uy')

                # changing the ticks on the y-axis
                ticks = fig.get_yticks() / _scale
                fig.set_yticklabels([round(x, 3) for x in ticks])

                _xdata = [p.x + dx * _scale for p, dx in zip(beam.nodes, dxs)]
                _ydata = [p.y + dy * _scale for p, dy in zip(beam.nodes, dys)]

                plt.plot(_xdata, _ydata, 'k-', linewidth=1, zorder=4)  # beams
                plt.scatter(_xdata, _ydata, marker='s', color='k', s=30, zorder=4)  # nodes

                for xy in zip(_xdata, _ydata):
                    fig.annotate('(%.3f, %.3f)' % (xy[0], xy[1] / _scale), xy=xy, textcoords='data', )

        if show:
            plt.axis('tight')
            plt.axis('equal')
            plt.show()
