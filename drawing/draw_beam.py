# -*- coding: utf-8 -*-
__author__ = 'jakabgabor'

from drawing import _plotting_available, plt


def draw_structure(structure, show=True, deformed=True):

    if _plotting_available:
        for beam in structure.beams:
            plt.plot([p.x for p in beam.nodes], [p.y for p in beam.nodes], 'k-')  # beams
            plt.plot([p.x for p in beam.nodes], [p.y for p in beam.nodes], 'r.', markersize=8)  # nodes

            # plot BCs (axis-style)
            # plot loads

            # http://12000.org/my_notes/stiffness_matrix/stiffness_matrix_report.htm


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
