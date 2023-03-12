import numpy as np
from math import sin, cos
import matplotlib.pyplot as plt


def RefineU(Ua, BC):
    '''
        Take active displacements and boundary conditions, return nodal displacements, {U}.

        Parameters
        ----------
        Ua : numpy array
            Active (Unconstrained) displacements matrix, {Ua}.
        BC : list
            Boundary conditions.

        Returns
        -------
        Unew : numpy array
            Nodal displacements, {U}.
    '''
    # nd = len(BC) * 2
    nd = len(BC)+Ua.shape[0]
    Unew = np.zeros(nd)
    BCnew = [(x[0]-1)*2+x[1]-1 for x in BC]
    BCval = [x[2] for x in BC]
    j,k = 0, 0
    for i in range(nd):
        if i in BCnew:
            Unew[i] = BCval[j]
            j += 1
        else:
            Unew[i] = Ua[k]
            k += 1
    return Unew


def RefineF(Fc, F):
    '''
        Take reaction forces and external forces, return nodal forces matrix, {F}.

        Parameters
        ----------
        Fc : numpy array
            Reaction forces of truss (constrained forces matrix, {Fc}).
        F : list
            External forces of truss.

        Returns
        -------
        Fnew : numpy array
            Nodal forces, {F}.
    '''
    nd = len(F) * 2
    Fnew = np.zeros(nd)
    Fnew = np.concatenate((Fc, Fnew))
    for i in F:
        Fnew[(i[0]*2)-2] = i[1]
        Fnew[(i[0]*2)-1] = i[2]
    return Fnew


def Strain(U, df):
    '''
        Take nodal displacements and dataframe of the problem, return elements' strain.

        Parameters
        ----------
        U : numpy array
            Nodal displacements.
        df : pandas dataframe
            Geometry definition of the problem.

        Returns
        -------
        df : pandas dataframe
            With strain of each element added to it.
    '''
    strain = []
    for i in df.index:
        node1 = df.iloc[i, 1][0]
        node2 = df.iloc[i, 1][1]
        Direction = df.iloc[i, 5]
        Length = df.iloc[i, 2]
        u1 = U[(2*node1)-2] * cos(Direction) + U[2*node1-1] * sin(Direction)
        u2 = U[(2*node2)-2] * cos(Direction) + U[2*node2-1] * sin(Direction)
        strain.append((u2-u1)/Length)
    df.insert(loc=6, column='Strain', value=strain)
    return df


def Stress(df):
    '''
        Take nodal displacements and dataframe of the problem, return elements' stress.

        Parameters
        ----------
        df : pandas dataframe
            Geometry definition of the problem.

        Returns
        -------
        df : pandas dataframe
            With stress of each element added to it.
    '''
    stress = []
    for i in df.index:
        E = df.iloc[i, 4]
        epsilon = df.iloc[i, 6]
        stress.append(E * epsilon)
    df.insert(loc=6, column='Stress', value=stress)
    return df


def ShowTruss(Coord, ElmCon):
    """
    Plot truss structure with node number and element number annotation.

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
    None.

    """

    x = [Coord[i][0] for i in range(len(Coord))]
    y = [Coord[i][1] for i in range(len(Coord))]

    ratio = (max(y)-min(y))/(max(x)-min(x))

    plt.figure(figsize=(10, 10*ratio), dpi=50)
    ax = plt.axes()
    plt.ylim(np.min(y) - np.max(y)*0.1, np.max(y) + np.max(y)*0.1)
    plt.xlim(np.min(x) - np.max(x)*0.1, np.max(x) + np.max(x)*0.1)
    ax.scatter(x, y)
    ma = np.max(x) * 0.01
    for num in range(len(Coord)):
        ax.annotate(str(num+1), (Coord[num][0] + ma, Coord[num][1] + ma),
                    bbox=dict(boxstyle="circle", alpha=0.1), fontsize=12)

    for i, sele in enumerate(ElmCon):
        ax.plot([Coord[sele[0]-1][0], Coord[sele[1]-1][0]],
                [Coord[sele[0]-1][1], Coord[sele[1]-1][1]], 'k', lw=1)
        ax.annotate(str(i+1), ((Coord[sele[0]-1][0]+Coord[sele[1]-1][0])/2 + ma, (Coord[sele[0]-1][1] +
                    Coord[sele[1]-1][1])/2 + ma), bbox=dict(boxstyle="square", fc="g", alpha=0.2), fontsize=12)


def ShowDeformedTruss(Coord, ElmCon, U, scale=10.0):
    """
    Plot truss structure with node number and element number annotation. 
    Also plots the deformed truss structure after loading.

    Parameters
    ----------
    Coord : list
        Nodal coordrindates. Each item is a tuple of (x,y)
        coordinates of nodes.
    ElmCon : list
        Element connectivity table. Each item is a tuple
        that contains first and second node number of bar.
    U : numpy array
        Nodal displacements.
    scale : float, optional
        Deformation scaling. The default is 10.0 .


    Returns
    -------
    None.

    """

    x = [Coord[i][0] for i in range(len(Coord))]
    y = [Coord[i][1] for i in range(len(Coord))]
    ratio = (max(y)-min(y))/(max(x)-min(x))

    plt.figure(figsize=(10, 10*ratio), dpi=50)
    ax = plt.axes()
    ax.scatter(x, y)
    ma = np.max(x) * 0.01
    for num in range(len(Coord)):
        ax.annotate(str(num+1), (x[num] + ma, y[num] + ma),
                    bbox=dict(boxstyle="circle", alpha=0.1), fontsize=12)

    for i, sele in enumerate(ElmCon):
        ax.plot([x[sele[0]-1], x[sele[1]-1]],
                [y[sele[0]-1], y[sele[1]-1]], 'k', lw=1)
        ax.annotate(str(i+1), ((x[sele[0]-1]+x[sele[1]-1])/2 + ma, (y[sele[0]-1]+y[sele[1]-1]
                                                                    )/2 + ma), bbox=dict(boxstyle="square", fc="g", alpha=0.2), fontsize=12)

    U = U.reshape((-1, 2))*scale
    newx = [tx+U[i, 0] for i, tx in enumerate(x)]
    newy = [ty+U[i, 1] for i, ty in enumerate(y)]
    ax.scatter(newx, newy)
    for i, sele in enumerate(ElmCon):
        ax.plot([newx[sele[0]-1], newx[sele[1]-1]],
                [newy[sele[0]-1], newy[sele[1]-1]], 'r:', lw=1)


def ShowTrussCross(Coord, ElmCon, A):
    """
    Plot truss structure with node number and element number annotation.
    Maps elements cross-sections to linewidth between [1,5]

    Parameters
    ----------
    Coord : list
        Nodal coordrindates. Each item is a tuple of (x,y)
        coordinates of nodes.
    ElmCon : list
        Element connectivity table. Each item is a tuple
        that contains first and second node number of bar.
    A : list
        Elements Cross-sections.

    Returns
    -------
    None.

    """

    mx, mn = max(A), min(A)
    b, a = 5, 1
    normalA = [(((b-a)*(x-mn))/(mx-mn))+a for x in A]

    x = [Coord[i][0] for i in range(len(Coord))]
    y = [Coord[i][1] for i in range(len(Coord))]

    ratio = (max(y)-min(y))/(max(x)-min(x))

    plt.figure(figsize=(10, 10*ratio), dpi=50)
    ax = plt.axes()
    plt.ylim(np.min(y) - np.max(y)*0.1, np.max(y) + np.max(y)*0.1)
    plt.xlim(np.min(x) - np.max(x)*0.1, np.max(x) + np.max(x)*0.1)
    ax.scatter(x, y)
    ma = np.max(x) * 0.01
    for num in range(len(Coord)):
        ax.annotate(str(num+1), (Coord[num][0] + ma, Coord[num][1] + ma),
                    bbox=dict(boxstyle="circle", alpha=0.1), fontsize=12)

    for i, sele in enumerate(ElmCon):
        ax.plot([Coord[sele[0]-1][0], Coord[sele[1]-1][0]],
                [Coord[sele[0]-1][1], Coord[sele[1]-1][1]], 'k', lw=normalA[i])
        ax.annotate(str(i+1), ((Coord[sele[0]-1][0]+Coord[sele[1]-1][0])/2 + ma, (Coord[sele[0]-1][1] +
                    Coord[sele[1]-1][1])/2 + ma), bbox=dict(boxstyle="square", fc="g", alpha=0.2), fontsize=12)


def PlotWDS(w, d, s):
    """
    Plot structure weight vs. maximum deflection vs. maximum stress where
    maximum stresses are mapped to markers size between [20,80].

    Parameters
    ----------
    w : list
        A list of weights.
    d : TYPE
        A list of maximum deflections.
    s : TYPE
        A list of maximum Stresses.

    Returns
    -------
    None.

    """
    plt.rcParams["figure.figsize"] = (5, 5)
    mx, mn = max(s), min(s)
    b, a = 80, 20
    nst = [(((b-a)*(x-mn))/(mx-mn))+a for x in s]
    plt.figure()
    ax = plt.axes()
    ax.scatter(w, d, s=nst, alpha=0.3)
    plt.xlabel('Weight (lb.)')
    plt.ylabel('Maximum Deflection (in.)')


def PlotWD(w, D):
    '''
    Plot structure weight vs. maximum deflection

    Parameters
    ----------
    w : list
        A list of weights.
    D : TYPE
        A list of maximum deflections.

    Returns
    -------
    None.

    '''
    plt.rcParams["figure.figsize"] = (8, 5*len(w))
    figure, axis = plt.subplots(len(w), 1)
    for j in range(len(w)):
        axis[j].plot(w[j], D[j])
        axis[j].set_title("Bar %s" % (j+1))
    for ax in axis.flat:
        ax.set(xlabel='Weight (lb.)', ylabel='Maximum Deflection (in.)')


def ShowTrussElast(Coord, ElmCon, E):
    """
    Plot truss structure with node number and element number annotation.
    Maps elements elasticity modulus to colors of black or red.

    Parameters
    ----------
    Coord : list
        Nodal coordrindates. Each item is a tuple of (x,y)
        coordinates of nodes.
    ElmCon : list
        Element connectivity table. Each item is a tuple
        that contains first and second node number of bar.
    E : list
        Elements' elastic modulus.

    Returns
    -------
    None.

    """

    mx, mn = max(E), min(E)
    b, a = 5, 1
    normalA = [(((b-a)*(x-mn))/(mx-mn))+a for x in E]
    for i in range(len(normalA)):
        if normalA[i] == 1:
            normalA[i] = 'r'
        else:
            normalA[i] = 'k'
    x = [Coord[i][0] for i in range(len(Coord))]
    y = [Coord[i][1] for i in range(len(Coord))]

    ratio = (max(y)-min(y))/(max(x)-min(x))

    plt.figure(figsize=(10, 10*ratio), dpi=50)
    ax = plt.axes()
    plt.ylim(np.min(y) - np.max(y)*0.1, np.max(y) + np.max(y)*0.1)
    plt.xlim(np.min(x) - np.max(x)*0.1, np.max(x) + np.max(x)*0.1)
    ax.scatter(x, y)
    ma = np.max(x) * 0.01
    for num in range(len(Coord)):
        ax.annotate(str(num+1), (Coord[num][0] + ma, Coord[num][1] + ma),
                    bbox=dict(boxstyle="circle", alpha=0.1), fontsize=12)

    for i, sele in enumerate(ElmCon):
        ax.plot([Coord[sele[0]-1][0], Coord[sele[1]-1][0]],
                [Coord[sele[0]-1][1], Coord[sele[1]-1][1]], normalA[i])
        ax.annotate(str(i+1), ((Coord[sele[0]-1][0]+Coord[sele[1]-1][0])/2 + ma, (Coord[sele[0]-1][1] +
                    Coord[sele[1]-1][1])/2 + ma), bbox=dict(boxstyle="square", fc="g", alpha=0.2), fontsize=12)
