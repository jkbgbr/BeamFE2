# -*- coding: utf-8 -*-


class Node(object):
    """
    Node objects to be used with the FE Model
    """

    def __init__(self, ID=None, coords=()):
        self.ID = ID
        self.coords = coords

    def __repr__(self):
        return 'Node(ID=%d, coords=(%.2f, %.2f)' % (self.ID, self.x, self.y)

    @classmethod
    def from_dict(cls, adict):
        try:
            node = Node(
                ID=adict['ID'],
                coords=adict['coords'],
            )
            assert all([hasattr(node, x) for x in adict.keys()])  # check if the dict is too long
        except KeyError as e:  # check if the dict is too short
            raise Exception('Missing key from dict when creating %s: %s' % (cls.__name__, e))
        return node

    def set_coords(self, coords):
        self.coords = coords

    @property
    def x(self):
        return self.coords[0]

    @property
    def y(self):
        return self.coords[1]
