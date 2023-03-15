# %% [0] Import libs:
import matplotlib.pyplot as plt
import myGenetic as GA
import numpy as np
import LoadCase1 as Inp
import Pre
import Sol
import Post
import time
# %% [1] Main module functions:
START = time.time()


def Run(Coord, ElmCon, A, E, theta, BC, F):
    """
    Perform finite element analysis of truss.

    Parameters
    ----------
    Coord : list
        Nodal coordrindates. Each item is a tuple of (x,y)
        coordinates of nodes.
    ElmCon : list
        Element connectivity table. Each item is a tuple
        that contains first and second node number of bar.
    A : list
        Element cross section area table. Each item represents the cross section area
        of bar. indices match with ElmCon.
    E : list
        Element elastic modulus table. Each item represents the elastic modulus
        of bar. indices match with ElmCon.
    theta : list
        Element direction table. Each item represents the angle between the bar and
        the global X axis in radians. indices match with ElmCon.
    BC : list
        Boundary conditions.
    F : list
        External forces of truss.

    Returns
    -------
    U : ndarray
        Nodal displacements, {U}.
    df : pandas dataframe
        Geometry definition of the problem.

    """

    L = Pre.LenCalc(Coord, ElmCon)
    df = Pre.dataframe(ElmCon, A, E, L, theta)
    df = Pre.ElmStf(df)
    K = Pre.TotStf(df)
    Ua = Sol.activeDisplacement(BC, K, F)
    # Fc = Sol.constrainedForces(BC, K, Ua)
    U = Post.RefineU(Ua, BC)
    # F = Post.RefineF(Fc, F)
    df = Post.Strain(U, df)
    df = Post.Stress(df)

    return U, df
# %% [2] Run an example:

def validation():
    
    Coord, ElmCon, A, E, BC, F = Inp.Coord, Inp.ElmCon, Inp.A, Inp.E, Inp.BC, Inp.F
    Post.ShowTruss(Coord, ElmCon)
    theta = Pre.angCalc(Coord, ElmCon)
    U, ansdf = Run(Coord, ElmCon, A, E, theta, BC, F)
    Post.ShowDeformedTruss(Coord, ElmCon, U)
    
    return U
# %% 10-bar planar truss:
Coord, ElmCon, E, BC, F = Inp.Coord, Inp.ElmCon, Inp.E, Inp.BC, Inp.F
theta = Pre.angCalc(Coord, ElmCon)
density, Length = Inp.density, Pre.LenCalc(Coord, ElmCon)

nd = 10
MaxDeflection = 2
MaxStress = 25000

A = GA.createPopReal(300, nd, [[10**-5, 35]])


def Weight(A, rho=density, L=Length):
    pen = 0
    U, ans = Run(Coord, ElmCon, A, E, theta, BC, F)
    if (max(abs(U)) > MaxDeflection):
        pen += (max(abs(U)) - MaxDeflection)/MaxDeflection
    if (max(abs(ans["Stress"])) > MaxStress):
        pen += (max(abs(ans["Stress"])) > MaxStress)/MaxStress

    f = sum([lA*lL*ldensity for lA, lL, ldensity in zip(A, L, rho)]) + pen*10**6

    return -f


nPar = 200
fitlist = []

for i in range(50):
    fitness = GA.fitCal(Weight, A)
    print(i)
    print("Best fitness = %f" % -max(fitness))
    fitlist.append(max(fitness))
    parents = GA.parSelRWS(A, fitness, nPar, rechoose=True)
    children = np.array(list(map(
        GA.combConWHA, parents[0:nPar//2, :], parents[nPar//2:nPar, :]))).reshape((nPar, nd))
    np.apply_along_axis(GA.mutConUni, 1, children, 0.09, 35.0, 0.2)
    A = GA.surSelRep(A, children, fitness, GA.fitCal(
        Weight, children), dup=True)


END = time.time()
print(END-START)
print(A[np.argsort(list(map(abs, GA.fitCal(Weight, A))))][:3])
print(np.sort(list(map(abs, GA.fitCal(Weight, A))))[0:3])


Post.ShowTrussCross(Coord, ElmCon, A[np.argsort(list(map(abs, GA.fitCal(Weight, A))))][0])
plt.figure()
plt.plot(fitlist)
