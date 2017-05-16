# -*- coding: utf-8 -*-

from Modell import HermitianBeam_2D as HB
from Modell import Structure
from Modell import Node
import math

"""
This example is from the opensees website: 
http://opensees.berkeley.edu/wiki/index.php/Eigen_analysis_of_a_two-story_shear_frame
origonal: Anil K. Chopra, Example in 10.5
the opensees tcl source is included at the end of this file
Kinda OK.
"""


h = 120  # story heigth
L = 2 * h  # bay width
mass_1 = 1000000 / 386.
rho = 0.01  # steel / 1000
EE = 29000
I = 320
A = 63.41

_nodes = [0]  # 0th element to have mathcing ID and list position
_nodes.append(Node.Node.from_dict(adict={'ID': 1, 'coords': (0, 0)}))
_nodes.append(Node.Node.from_dict(adict={'ID': 2, 'coords': (L, 0)}))
_nodes.append(Node.Node.from_dict(adict={'ID': 3, 'coords': (0, h)}))
_nodes.append(Node.Node.from_dict(adict={'ID': 4, 'coords': (L, h)}))
_nodes.append(Node.Node.from_dict(adict={'ID': 5, 'coords': (0, 2 * h)}))
_nodes.append(Node.Node.from_dict(adict={'ID': 6, 'coords': (L, 2 * h)}))

# beams
_beams = []
_beams.append(HB.HermitianBeam2D.from_dict(adict={'ID': 1, 'E': EE, 'I': 2 * I, 'A': 1 * A, 'rho': rho, 'i': _nodes[1], 'j': _nodes[3]}))
_beams.append(HB.HermitianBeam2D.from_dict(adict={'ID': 2, 'E': EE, 'I': 1 * I, 'A': 1 * A, 'rho': rho, 'i': _nodes[3], 'j': _nodes[5]}))

_beams.append(HB.HermitianBeam2D.from_dict(adict={'ID': 3, 'E': EE, 'I': 2 * I, 'A': 1 * A, 'rho': rho, 'i': _nodes[2], 'j': _nodes[4]}))
_beams.append(HB.HermitianBeam2D.from_dict(adict={'ID': 4, 'E': EE, 'I': 1 * I, 'A': 1 * A, 'rho': rho, 'i': _nodes[4], 'j': _nodes[6]}))

_beams.append(HB.HermitianBeam2D.from_dict(adict={'ID': 5, 'E': EE, 'I': 2 * I, 'A': I, 'rho': rho, 'i': _nodes[3], 'j': _nodes[4]}))
_beams.append(HB.HermitianBeam2D.from_dict(adict={'ID': 6, 'E': EE, 'I': 2 * I, 'A': I, 'rho': rho, 'i': _nodes[5], 'j': _nodes[6]}))

# supports
BCs = {1: ['ux', 'uy', 'rotz'], 2: ['ux', 'uy', 'rotz']}  # supports as dict

# this is the cantilever itself, composed of the beams, complete with supports
structure = Structure.Structure(beams=_beams, supports=BCs)

# adding masses for the moda lanalysis
structure.add_nodal_mass(nodeID=3, mass={'mx': mass_1, 'my': mass_1})
structure.add_nodal_mass(nodeID=4, mass={'mx': mass_1, 'my': mass_1})
structure.add_nodal_mass(nodeID=5, mass={'mx': mass_1/2., 'my': mass_1/2.})
structure.add_nodal_mass(nodeID=6, mass={'mx': mass_1/2., 'my': mass_1/2.})

# adding loads; these are not taken into account when performing the modal analysis
for b in structure.beams:
    structure.add_internal_loads(beam=b, loadtype='uniform perpendicular force', value=5.00)

# solving it
structure.solver['linear static'].solve()
structure.draw(analysistype='linear static')

structure.solver['modal'].solve()

# posprocessing
om_1 = (2.198 * math.sqrt(EE * I / (mass_1 * h ** 3))) / (2 * math.pi)
om_2 = (5.85 * math.sqrt(EE * I / (mass_1 * h ** 3))) / (2 * math.pi)
print('expected / calculated values in [Hz]:')
print('Mode 1: %.8f / %.8f' % (om_1, structure.results['modal'].frequencies[0]))
print('Mode 2: %.8f / %.8f' % (om_2, structure.results['modal'].frequencies[1]))

for i in range(2):
    structure.draw(analysistype='modal', mode=i)


# # Eigen analysis of a two-storey one-bay frame; Example 10.5 from "Dynamics of Structures" book by Anil Chopra
#
# # units: kips, in, sec
#
# #       Vesna Terzic, 2010
#
# #delete all previosly constructed objects
# wipe;
#
# #set input variables
# #--------------------
#
# #mass
# set m  [expr 100.0/386.0]
#
# #number of modes
# set numModes 2
#
# #material
# set A 63.41
# set I 320.0
# set E 29000.0
#
# #geometry
# set L 240.
# set h  120.
#
# # create data directory
# file mkdir modes;
#
# # define the model
# #---------------------------------
# #model builder
# model BasicBuilder -ndm 2 -ndf 3
#
# # nodal coordinates:
# node 1   0.  0. ;
# node 2   $L  0. ;
# node 3   0.  $h ;
# node 4   $L  $h ;
# node 5   0.  [expr 2*$h];
# node 6   $L  [expr 2*$h];
#
# # Single point constraints -- Boundary Conditions
# fix 1 1 1 1;
# fix 2 1 1 1;
#
# # assign mass
# mass 3 $m 0. 0. ;
# mass 4 $m 0. 0. ;
# mass 5 [expr $m/2.] 0. 0. ;
# mass 6 [expr $m/2.] 0. 0. ;
#
# # define geometric transformation:
# set TransfTag 1;
# geomTransf Linear $TransfTag ;
#
# # define elements:
# # columns
# element elasticBeamColumn 1 1 3 $A $E [expr 2.*$I] $TransfTag;
# element elasticBeamColumn 2 3 5 $A $E $I           $TransfTag;
# element elasticBeamColumn 3 2 4 $A $E [expr 2.*$I] $TransfTag;
# element elasticBeamColumn 4 4 6 $A $E $I           $TransfTag;
# # beams
# element elasticBeamColumn 5 3 4 $A $E [expr 2.*$I] $TransfTag;
# element elasticBeamColumn 6 5 6 $A $E $I           $TransfTag;
#
# # record eigenvectors
# #----------------------
# for { set k 1 } { $k <= $numModes } { incr k } {
#     recorder Node -file [format "modes/mode%i.out" $k] -nodeRange 1 6 -dof 1 2 3  "eigen $k"
# }
#
# # perform eigen analysis
# #-----------------------------
# set lambda [eigen  $numModes];
#
# # calculate frequencies and periods of the structure
# #---------------------------------------------------
# set omega {}
# set f {}
# set T {}
# set pi 3.141593
#
# foreach lam $lambda {
# 	lappend omega [expr sqrt($lam)]
# 	lappend f [expr sqrt($lam)/(2*$pi)]
# 	lappend T [expr (2*$pi)/sqrt($lam)]
# }
#
# puts "periods are $T"
#
# # write the output file cosisting of periods
# #--------------------------------------------
# set period "modes/Periods.txt"
# set Periods [open $period "w"]
# foreach t $T {
# 	puts $Periods " $t"
# }
# close $Periods
#
# # record the eigenvectors
# #------------------------
# record
#
# # create display  for mode shapes
# #---------------------------------
# #                 $windowTitle $xLoc $yLoc $xPixels $yPixels
# recorder display "Mode Shape 1"  10    10     500      500     -wipe
# prp $h $h 1;                                         # projection reference point (prp); defines the center of projection (viewer eye)
# vup  0  1 0;                                         # view-up vector (vup)
# vpn  0  0 1;                                         # view-plane normal (vpn)
# viewWindow -200 200 -200 200;                        # coordiantes of the window relative to prp
# display -1 5 20;                                     # the 1st arg. is the tag for display mode (ex. -1 is for the first mode shape)
#                                                      # the 2nd arg. is magnification factor for nodes, the 3rd arg. is magnif. factor of deformed shape
# recorder display "Mode Shape 2" 10 510 500 500 -wipe
# prp $h $h 1;
# vup  0  1 0;
# vpn  0  0 1;
# viewWindow -200 200 -200 200
# display -2 5 20
#
# # get values of eigenvectors for translational DOFs
# #---------------------------------------------------
# set f11 [nodeEigenvector 3 1 1]
# set f21 [nodeEigenvector 5 1 1]
# set f12 [nodeEigenvector 3 2 1]
# set f22 [nodeEigenvector 5 2 1]
# puts "eigenvector 1: [list [expr {$f11/$f21}] [expr {$f21/$f21}] ]"
# puts "eigenvector 2: [list [expr {$f12/$f22}] [expr {$f22/$f22}] ]"
