# -*- coding: utf-8 -*-

from drawing import _plotting_available, plt
from Modell import HermitianBeam_2D as HeBe


def draw_structure(structure, show=True, analysistype=None, mode=0):
    if _plotting_available:
        for beam in structure.beams:
            # line width is in pixel
            plt.plot([p.x for p in beam.nodes], [p.y for p in beam.nodes], 'royalblue', linewidth=3, alpha=1, zorder=1)  # beams
            plt.scatter([p.x for p in beam.nodes], [p.y for p in beam.nodes], marker='s', color='navy', alpha=0.5, s=30, zorder=2)  # nodes

        # plot supports
        # 1/scale will be the length of the bars representing the restricted translational degrees of freedom
        scale = 20
        figsize = [plt.rcParams["figure.dpi"] * x for x in plt.rcParams["figure.figsize"]]  # width of the fig in pixels
        w_figsize = figsize[0]
        supportsize = w_figsize / scale
        for k, v in structure.supports.items():
            _aktnode = [x for x in structure.nodes if x.ID == k][0]
            for dof in v:
                if dof == 'ux':  # horizonatal
                    plt.plot([_aktnode.x, _aktnode.x + supportsize], [_aktnode.y, _aktnode.y], 'g-', linewidth=4, zorder=6)  # a horizontal line
                if dof == 'uy':  # vertical
                    plt.plot([_aktnode.x, _aktnode.x], [_aktnode.y, _aktnode.y - supportsize], 'g-', linewidth=4, zorder=6)  # a horizontal line
                if dof == 'rotz':  # rotation about Z
                    plt.plot([_aktnode.x, _aktnode.x], [_aktnode.y, _aktnode.y], 'seagreen', markersize=16, zorder=7)  # a point

        # plot the deformed shape
        # length of the longes beam - this will be the base for the scaling
        _long = sorted(structure.beams, key=lambda x: x.l)[-1].l

        # # displacements

        dre = structure.results[analysistype].global_displacements(mode=0, asvector=True)
        # # the scaling factor, based on the larges displacement and the length of the longest element

        _scale = _long / max(abs(dre))

        for beam in structure.beams:
            # beam displacements by component
            dxs = structure.results[analysistype].element_displacements(mode=mode, beam=beam)['ux']
            dys = structure.results[analysistype].element_displacements(mode=mode, beam=beam)['uy']

            # data to plot: original positions + displacements
            _xdata = [p.x + dx * _scale for p, dx in zip(beam.nodes, dxs)]
            _ydata = [p.y + dy * _scale for p, dy in zip(beam.nodes, dys)]

            # the nodes as squares
            plt.scatter(_xdata, _ydata, marker='s', color='k', s=30, zorder=3)

            # plot the deformed shape - using the internal points from the shape functions
            _deflected = beam.deflected_shape(local=False, scale=_scale, disps=structure.results[analysistype].element_displacements(mode=mode, beam=beam, asvector=True))
            plt.plot([x[0] for x in _deflected], [x[1] for x in _deflected], 'k-', zorder=3)

        # # plot loads - concentrated forces only for now
        # for lindex, load in enumerate(HeBe.np_matrix_tolist(structure.q)):
        #     if load != 0:
        #         _node, component = structure.nodenr_dof_from_position(position=lindex)
        #         mp = structure.node_by_ID(id=_node).coords  # starting point of the arrow
        #
        #         if component == 'FX':
        #             _norm = (load * supportsize / abs(load), 0,)
        #         elif component == 'FY':
        #             _norm = (0, load * supportsize / abs(load))
        #         else:
        #             _norm = False
        #         # plotting, if there is a norm
        #         ax = plt.gca()
        #         if _norm:
        #             ax.arrow(mp[0], mp[1], _norm[0], _norm[1], head_width=supportsize/20., head_length=supportsize/20., fc='r', ec='r')

        if show:
            plt.axis('tight')
            plt.axis('equal')
            plt.show()

        return True

    else:
        return False


def _draw_structure(structure, show=True, displacements=True, internal_actions=False):

    if _plotting_available:
        for beam in structure.beams:
            # line width is in pixel
            plt.plot([p.x for p in beam.nodes], [p.y for p in beam.nodes], 'royalblue', linewidth=3, alpha=1, zorder=1)  # beams
            plt.scatter([p.x for p in beam.nodes], [p.y for p in beam.nodes], marker='s', color='navy', alpha=0.5, s=30, zorder=2)  # nodes

        # plot supports
        # 1/scale will be the length of the bars representing the restricted translational degrees of freedom
        scale = 20
        figsize = [plt.rcParams["figure.dpi"] * x for x in plt.rcParams["figure.figsize"]]  # width of the fig in pixels
        w_figsize = figsize[0]
        supportsize = w_figsize / scale
        for k, v in structure.supports.items():
            _aktnode = [x for x in structure.nodes if x.ID == k][0]
            for dof in v:
                if dof == 'ux':  # horizonatal
                    plt.plot([_aktnode.x, _aktnode.x + supportsize], [_aktnode.y, _aktnode.y], 'g-', linewidth=4, zorder=6)  # a horizontal line
                if dof == 'uy':  # vertical
                    plt.plot([_aktnode.x, _aktnode.x], [_aktnode.y, _aktnode.y - supportsize], 'g-', linewidth=4, zorder=6)  # a horizontal line
                if dof == 'rotz':  # rotation about Z
                    plt.plot([_aktnode.x, _aktnode.x], [_aktnode.y, _aktnode.y], 'seagreen', markersize=16, zorder=7)  # a point

        # plot the deformed shape
        if displacements:
            ax = plt.gca()
            # length of the longes beam - this will be the base for the scaling
            _long = sorted(structure.beams, key=lambda x: x.l)[-1].l

            # displacements
            dre = structure.resulting_displacement
            # the scaling factor, based on the larges displacement and the length of the longest element
            _scale = (_long / 10.) / max(dre)

            for beam in structure.beams:
                # beam displacements by component
                dxs = beam.displacement_component(component='ux')
                dys = beam.displacement_component(component='uy')

                # data to plot: original positions + displacements
                _xdata = [p.x + dx * _scale for p, dx in zip(beam.nodes, dxs)]
                _ydata = [p.y + dy * _scale for p, dy in zip(beam.nodes, dys)]

                # the nodes as squares
                plt.scatter(_xdata, _ydata, marker='s', color='k', s=30, zorder=3)

                # plot the deformed shape - using the internal points from the shape functions
                _deflected = beam.deflected_shape(local=False, scale=_scale)
                plt.plot([x[0] for x in _deflected], [x[1] for x in _deflected], 'k-', zorder=3)

        if internal_actions:
            pass

        # plot loads - concentrated forces only for now
        for lindex, load in enumerate(HeBe.np_matrix_tolist(structure.q)):
            if load != 0:
                _node, component = structure.nodenr_dof_from_position(position=lindex)
                mp = structure.node_by_ID(id=_node).coords  # starting point of the arrow

                if component == 'FX':
                    _norm = (load * supportsize / abs(load), 0,)
                elif component == 'FY':
                    _norm = (0, load * supportsize / abs(load))
                else:
                    _norm = False
                # plotting, if there is a norm
                ax = plt.gca()
                if _norm:
                    ax.arrow(mp[0], mp[1], _norm[0], _norm[1], head_width=supportsize/20., head_length=supportsize/20., fc='k', ec='k')

        if show:
            plt.axis('tight')
            plt.axis('equal')
            plt.show()

        return True

    else:
        return False
