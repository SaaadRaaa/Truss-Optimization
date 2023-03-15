from math import pi

Coord = [(0, 0), (0, 360), (360, 0), (360, 360), (720, 0), (720, 360)]

ElmCon = [(1, 3), (1, 4), (2, 3), (2, 4), (3, 4),
          (3, 5), (3, 6), (4, 5), (4, 6), (5, 6)]

crsSctn = 25
A = [crsSctn, crsSctn, crsSctn, crsSctn, crsSctn,
     crsSctn, crsSctn, crsSctn, crsSctn, crsSctn]

ElasMod = 10**7
E = [ElasMod, ElasMod, ElasMod, ElasMod, ElasMod,
     ElasMod, ElasMod, ElasMod, ElasMod, ElasMod]

# [number of node, DoF(1 or 2), value]
BC = [[1, 1, 0], [1, 2, 0], [2, 1, 0], [2, 2, 0]]

P1 = 150 * 1000
P2 = 50 * 1000
F = [[3, 0, -P1], [4, 0, P2], [5, 0, -P1], [6, 0, P2]]

rho = 0.1
density = [rho, rho, rho, rho, rho, rho, rho, rho, rho, rho]
