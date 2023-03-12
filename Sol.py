import numpy as np


def activeDisplacement(BC, K, F):
    """
        Take Boundary conditions, global stiffness matrix and forces. return active displacements, {Ua}.

        Parameters
        ----------
        BC : list
            Boundary conditions.
        K : ndarray
            Global stiffness matrix.
        F : list
            External forces of truss.

        Returns
        -------
        Ua : ndarray
            Active (Unconstrained) displacements matrix, {Ua}.
    """

    BCnew = [(x[0]-1)*2+x[1]-1 for x in BC]
    BCval = [x[2] for x in BC]
    Fnew = np.zeros(K.shape[0])
    
    for i,bc in enumerate(BCnew):
        Fnew += BCval[i]*K[:,bc]
    for node in F:
        Fnew[(node[0]-1)*2] += node[1] 
        Fnew[(node[0]-1)*2+1] += node[2] 
    Ft = np.delete(Fnew, BCnew)
    Knew = np.delete(K, BCnew, axis=0)
    Knew = np.delete(Knew, BCnew, axis=1)


    Ua = np.linalg.solve(Knew, Ft)
    return Ua

#%% Need to change 
#
#
# def constrainedForces(BC, K, Ua):
#     '''
#         Take Boundary conditions, global stiffness matrix and active displacements. return reaction forces, {Fc}.

#         Parameters
#         ----------
#         BC : list
#             Boundary conditions.
#         K : numpy array
#             Global stiffness matrix.
#         Ua : numpy array
#             Active (Unconstrained) displacements matrix, {Ua}.

#         Returns
#         -------
#         Fc : numpy array
#             Reaction forces of truss (constrained forces matrix, {Fc}).
#     '''
#     nd = K.shape[0]
#     dele = []
#     # BCnew1 = [x[0]*2-1 for x in BC]
#     # BCnew2 = [x[0]*2-2 for x in BC]
#     # BCnew = BCnew1 + BCnew2
#     BCnew = [(x[0]-1)*2+x[1]-1 for x in BC]
#     for x in range(nd):
#         if not x in BCnew:
#             dele.append(x)
#     Knew = np.delete(K, BCnew, axis=1)
#     Knew = np.delete(Knew, dele, axis=0)
#     Fc = np.matmul(Knew, Ua)
#     return Fc
#%%