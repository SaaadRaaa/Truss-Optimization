import numpy as np
import pandas as pd
from math import sin, cos, pi, atan


def LenCalc(Coord, ElmCon):
    """
    Calculate the lengths of bar elements.

    Parameters
    ----------
    Coord : list
        Nodal coordrindates. Each item is a tuple of (x,y)
        coordinates of nodes.
    ElmCon : list
        Element connectivity table. Each item is a tuple
        that contains first and second node number of bar.

    Returns
    -------
    L : list
        Element length table. Each item represents the length
        of bar. indices match with ElmC.
    """
    L = list()
    for elem in ElmCon:
        slength = (Coord[elem[0]-1][0]-Coord[elem[1]-1][0])**2 +\
                  (Coord[elem[0]-1][1]-Coord[elem[1]-1][1])**2
        L.append(np.sqrt(slength))

    return L


def angCalc(Coord, ElmCon):
    """
    Calculate the orientations of bar elements.

    Parameters
    ----------
    Coord : list
        Nodal coordrindates. Each item is a tuple of (x,y)
        coordinates of nodes.
    ElmCon : list
        Element connectivity table. Each item is a tuple
        that contains first and second node number of bar.

    Returns
    -------
    angs : list
        Element direction table. Each item represents the angle between the bar and
        the global X axis in radians. indices match with ElmCon.

    """

    angs = list()
    for elem in ElmCon:
        n1 = Coord[elem[0]-1]
        n2 = Coord[elem[1]-1]
        if n1[0] == n2[0]:
            angs.append(pi/2)
        elif n1[1] == n2[1]:
            angs.append(0)
        else:
            angs.append(atan((n2[1]-n1[1])/(n2[0]-n1[0])))

    return angs


def dataframe(ElmC, A, E, L, theta):
    """
    Take each element's connectivity, cross section area, elastic modulus, length and direction,
    return a pandas dataframe containing geometry definition of the problem.

    Parameters
    ----------
    ElmC : list
        Element connectivity table. Each item is a tuple
        that contains first and second node number of bar.
    A : list
        Element cross section area table. Each item represents the cross section area
        of bar. indices match with ElmC.
    E : list
        Element elastic modulus table. Each item represents the elastic modulus
        of bar. indices match with ElmC.
    L : list
        Element length table. Each item represents the length
        of bar. indices match with ElmC.
    theta : list
        Element direction table. Each item represents the angle between the bar and
        the global X axis in radians. indices match with ElmC.

    Returns
    -------
    df : pandas dataframe
        Geometry definition of the problem.
    """

    d = {'#': list(range(1, len(ElmC)+1)), 'Connectivity': ElmC, 'Length': L, 'Cross section area': A, 'Young\'s modulus': E,
         'Direction (Rad)': theta}
    df = pd.DataFrame(data=d)
    return df


def ElmStf(df):
    """
    Take problem's dataframe and add a column containing each element's
    stiffness matrix.

    Parameters
    ----------
    df : pandas dataframe
        Geometry definition of the problem.

    Returns
    -------
    df : pandas dataframe
        With element stiffness matrices added to it.
    """
    ke = []
    for i in df.index:
        k = df.iloc[i, 4] * df.iloc[i, 3] / df.iloc[i, 2]
        ke.append(np.zeros((4, 4)))
        ke[i][0, 0] = ke[i][2, 2] = k * cos(df.iloc[i, 5]) ** 2
        ke[i][1, 1] = ke[i][3, 3] = k * sin(df.iloc[i, 5]) ** 2
        ke[i][0, 1] = ke[i][1, 0] = k * cos(df.iloc[i, 5]) * sin(df.iloc[i, 5])
        ke[i][2, 3] = ke[i][3, 2] = k * cos(df.iloc[i, 5]) * sin(df.iloc[i, 5])
        ke[i][0, 2] = ke[i][2, 0] = -k * cos(df.iloc[i, 5]) ** 2
        ke[i][1, 3] = ke[i][3, 1] = -k * sin(df.iloc[i, 5]) ** 2
        ke[i][3, 0] = ke[i][2, 1] = -k * \
            cos(df.iloc[i, 5]) * sin(df.iloc[i, 5])
        ke[i][1, 2] = ke[i][0, 3] = -k * \
            cos(df.iloc[i, 5]) * sin(df.iloc[i, 5])
    df.insert(loc=6, column='Stiffness matrix', value=ke)
    return df


def TotStf(df):
    """
    Take problem's dataframe and return global stiffness matrix.

    Parameters
    ----------
    df : pandas dataframe
        Geometry definition of the problem.

    Returns
    -------
    K : ndarray
        Global stiffness matrix.
    """
    nd = max(df['Connectivity'].max())
    K = np.zeros((2*nd, 2*nd))
    for num in df.index:
        i = df.iloc[num, 1][0]-1
        j = df.iloc[num, 1][1]-1
        K[2*i, 2*i] += df.iloc[num, 6][0, 0]
        K[2*i+1, 2*i+1] += df.iloc[num, 6][1, 1]
        K[2*i+1, 2*i] += df.iloc[num, 6][0, 1]
        K[2*i, 2*i+1] += df.iloc[num, 6][1, 0]
        K[2*j, 2*j] += df.iloc[num, 6][2, 2]
        K[2*j+1, 2*j+1] += df.iloc[num, 6][3, 3]
        K[2*j+1, 2*j] += df.iloc[num, 6][2, 3]
        K[2*j, 2*j+1] += df.iloc[num, 6][3, 2]
        K[2*j, 2*i] += df.iloc[num, 6][0, 2]
        K[2*j+1, 2*i+1] += df.iloc[num, 6][1, 3]
        K[2*j+1, 2*i] += df.iloc[num, 6][3, 0]
        K[2*j, 2*i+1] += df.iloc[num, 6][2, 1]
        K[2*i, 2*j] += df.iloc[num, 6][2, 0]
        K[2*i+1, 2*j+1] += df.iloc[num, 6][3, 1]
        K[2*i+1, 2*j] += df.iloc[num, 6][1, 2]
        K[2*i, 2*j+1] += df.iloc[num, 6][0, 3]
    return K


def TrussWeight(A, L, density):
    """
    Calculate truss weight.

    Parameters
    ----------
    A : list
        A list of elements cross-sections.
    L : list
        A list of elements weights.
    density : list
        Material density.

    Returns
    -------
    _ : float
        Structure Weight.
    """
    return sum([lA*lL*ldensity for lA, lL, ldensity in zip(A, L, density)])
