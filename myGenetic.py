import numpy as np


# -----------------------------------------
# -----------------------------------------
# Creating Initial Population
# -----------------------------------------
# -----------------------------------------


def createPopReal(popSize, nDim=1, limCon=[[0, 4]]):
    """
    Create initial population with real value genes.

    Parameters
    ----------
    popSize : int
        Population size.
    nDim : int, default = 1
        Number of dimensions of the problem.
    limCon : list, default = [[0, 4]]
        Lower(inclusive) and upper(exclusive) bounds for each gene.
    """
    if len(limCon) < nDim:
        while len(limCon) < nDim:
            limCon.append(limCon[-1])

    rng = np.random.default_rng()
    initPop = rng.random((popSize, nDim))
    for i in range(nDim):
        a = limCon[i][0]
        b = limCon[i][1]
        initPop[:, i] = (b-a) * initPop[:, i] + a

    return initPop


def createPopInt(popSize, nDim=1, limCon=[[0, 4]]):
    """
    Create initial population with integer value genes.

    Parameters
    ----------
    popSize : int
        Population size.
    nDim : int, default = 1
        Number of dimensions of problem.
    limCon : list, default = [[0,4]] 
        Lower(inclusive) and upper(exclusive) bounds for each gene.
    """
    if len(limCon) < nDim:
        while len(limCon) < nDim:
            limCon.append(limCon[-1])

    rng = np.random.default_rng()
    lb = [con[0] for con in limCon]
    ub = [con[1] for con in limCon]
    initPop = rng.integers(lb, ub, size=(popSize, len(lb)))

    return initPop


def createPopPerm(popSize, ub, lb=1):
    """
    Create initial population with permutated value genes.

    Parameters
    ----------
    popSize : int
        Population size.
    ub : int
        Upper bound for values.
    lb : int, default = 1
        Lower bound for values.
    """
    rng = np.random.default_rng()
    arr = np.arange(lb, ub)
    initPop = rng.permutation(arr)
    for i in range(popSize-1):
        initPop = np.vstack([initPop, rng.permutation(arr)])

    return initPop


def createPopBin(popSize, nDim=4):
    """
    Create initial population with binary value genes.

    Parameters
    ----------
    popSize : int
        Population size.
    nDim : int, default = 4
        Number of dimensions of problem.
    """

    rng = np.random.default_rng()

    initPop = rng.integers(2, size=(popSize, nDim))

    return initPop


# -----------------------------------------
# -----------------------------------------
# Fitness Calculation
# -----------------------------------------
# -----------------------------------------

def fitCal(func, x):
    """
    Calculate fitness of one or more chromosomes.

    Parameters
    ----------
    func : function
        Fitness function.
    x : numpy array
        An array that each row represents a chromosome.

    Returns
    ----------
    _ : list
        A list of calculated fitnesses.
    """
    return [func(row) for row in x]


# -----------------------------------------
# -----------------------------------------
# Parent Selection
# -----------------------------------------
# -----------------------------------------


def fitScale(fitness, a=None, b=None):
    """
    Scale the fitness array with the expression a*fitness(i)+b and 
    return the normalized fitness. If the factors (a and b) are not given, 
    ensure that there would be no negative fitness value.

    Parameters
    ----------
    fitness : list
        A list containing parents fitness.
    a : float, default = None
    b : float, default = None
    """
    if a == None:
        a = 1.0
    if b == None:
        b = min(fitness)
        if b < 0:
            b = abs(b)
        else:
            b = 0

    newFitness = [a*x+b for x in fitness]
    sumFit = sum(newFitness)
    fitness = [x/sumFit for x in newFitness]

    return fitness


def parSelUNI(totPop, nPar=1, rechoose=False):
    """
        Random parent selection.

        Parameters
        ----------
        totPop : numpy array
            An array that each row represent a candidate parent.
        nPar : int, default = 1
            The number of desired parents.
        rechoose : boolean, default = False
            Wether the parent have the chance to be selected twice or not.

    """
    rng = np.random.default_rng()

    if rechoose == True:
        parInd = rng.integers(0, totPop.shape[0], size=1)
        for i in range(nPar-1):
            parInd = np.hstack(
                (parInd, rng.integers(0, totPop.shape[0], size=1)))
    else:
        parInd = rng.integers(0, totPop.shape[0], size=nPar)

    return totPop[parInd, :]


def parSelRWS(totPop, fitness, nPar=1, rechoose=False):
    """
    Fitness proportional parent selection.

    Parameters
    ----------
    totPop : numpy array
        An array that each row represent a candidate parent.
    fitness : list
        A list containing parents fitness.
    nPar : int, default = 1
        The number of desired parents.
    rechoose : boolean, default = False
        Wether the parent have the chance to be selected twice or not.

    """

    fitness = fitScale(fitness)
    totFit = sum(fitness)
    prob = [x/totFit for x in fitness]
    CDF = [sum(prob[0:i]) for i in range(1, len(prob)+1)]
    rng = np.random.default_rng()
    parInd = list()
    while len(parInd) < nPar:
        crus = rng.random()
        for i in range(len(CDF)):
            if crus > CDF[i]:
                continue
            else:
                ind = i
                break
        if rechoose == True:
            parInd.append(ind)
        elif ind not in parInd:
            parInd.append(ind)

    return totPop[parInd, :]


def parSelRNK(totPop, fitness, nPar=1, rechoose=False, probFunc='lin'):
    """
    Fitness ranking parent selection. This function use 'parSelRWS' function.

    Parameters
    ----------
    totPop : numpy array
        An array that each row represent a candidate parent.
    fitness : list
        A list containing parents fitness.
    nPar : int, default = 1
        The number of desired parents.
    rechoose : boolean, default = False
        Wether the parent have the chance to be selected twice or not.
    probFunc : 'lin' or 'exp', default = 'lin'
        Probability function.

    """

    fitness = fitScale(fitness)
    fitCopy = fitness.copy()
    inds = [fitCopy.index(x) for x in sorted(fitness)]
    if probFunc == 'lin':
        for i, ind in enumerate(inds):
            fitness[ind] = i+1
    elif probFunc == 'exp':
        for i, ind in enumerate(inds):
            fitness[ind] = np.exp(i+1)

    return parSelRWS(totPop, fitness, nPar, rechoose)


def parSelTNS(totPop, fitness, nPar=1, rechoose=False, indivs=3):
    """
    Hunger games!

    Parameters
    ----------
    totPop : numpy array
        An array that each row represent a candidate parent.
    fitness : list
        A list containing parents fitness.
    nPar : int, default = 1
        The number of desired parents.
    rechoose : boolean, default = False
        Wether the parent have the chance to be selected twice or not.
    indivs : int, default = 3
        The number of individuals participating in the tournaments.
    """
    if totPop.shape[0] <= 3:
        exit()

    selected_parents = []
    n = len(selected_parents)
    while n < nPar:
        rng = np.random.default_rng()
        i = rng.integers(0, totPop.shape[0], size=indivs)
        selected_indivs = totPop[i]
        selected_fits = np.array(fitness)[i]
        sorted_idx = np.argsort(-selected_fits)
        winner = selected_indivs[sorted_idx[0]].tolist()
        if rechoose == True:
            selected_parents.append(winner)
        elif rechoose == False and winner not in selected_parents:
            selected_parents.append(winner)
        n = len(selected_parents)
    parents = np.array([np.array(xi) for xi in selected_parents])

    return parents


# -----------------------------------------
# -----------------------------------------
# CrossOvers
# -----------------------------------------
# -----------------------------------------


def combDisOPC(P1, P2):
    """
    One point crossover for discrete values.
    """

    rng = np.random.default_rng()
    cInd = rng.integers(1, P1.shape[0], size=1)
    for i in range(int(cInd)):
        P1[i], P2[i] = P2[i], P1[i]

    return P1, P2


def combDisNPC(P1, P2, n=2):
    """
    N-point crossover for discrete values.
    """

    rng = np.random.default_rng()
    cInd = rng.integers(1, P1.shape[0], size=1)
    cInd2 = rng.integers(cInd, P1.shape[0], size=1)
    for i in range(int(cInd), int(cInd2)):
        P1[i], P2[i] = P2[i], P1[i]

    return P1, P2


def combDisUNC(P1, P2):
    """
    Uniform crossover for discrete values.
    """

    rng = np.random.default_rng()
    prob = rng.random(size=(1, P1.shape[0]))
    crit = rng.random()
    for i in range(P1.shape[0]):
        if prob[0, i] > crit:
            P1[i], P2[i] = P2[i], P1[i]

    return P1, P2


# ----------------------------------------


def combConSPA(P1, P2, alpha=0.5):
    """
    Simple Arithmetic combination for continuous values.
    """
    rng = np.random.default_rng()
    cInd = rng.integers(1, P1.shape[0], size=1)
    for i in range(int(cInd)):
        x = alpha*P1[i] + (1-alpha)*P2[i]
        P1[i], P2[i] = x, x

    return P1, P2


def combConSGA(P1, P2, alpha=0.5):
    """
    Single Arithmetic combination for continuous values.
    """

    rng = np.random.default_rng()
    cInd = rng.integers(0, P1.shape[0]+1, size=1)

    x = alpha*P1[cInd] + (1-alpha)*P2[cInd]
    P1[cInd], P2[cInd] = x, x

    return P1, P2


def combConWHA(P1, P2, alpha=0.5):
    """
    Whole arithmetic combination for continuous values.
    """

    ch1 = P1.copy()
    ch2 = P2.copy()

    for i in range(P1.shape[0]):
        ch1[i] = alpha*P1[i] + (1-alpha)*P2[i]
        ch2[i] = alpha*P2[i] + (1-alpha)*P1[i]

    return ch1, ch2


# ----------------------------------------


def combPerPMX(P1, P2):
    """
    Partially mapped crossover for permutated values.
    """

    ch1 = np.zeros_like(P1)
    ch2 = np.zeros_like(P1)

    rng = np.random.default_rng()
    cInd = rng.integers(1, P1.shape[0]+1, size=1)
    cInd2 = rng.integers(cInd, P1.shape[0], size=1)
    for i in range(int(cInd), int(cInd2)):
        ch1[i], ch2[i] = P1[i], P2[i]

    for i in range(int(cInd), int(cInd2)):
        if P2[i] not in ch1:
            tInd = i
            while True:
                tInd = np.where(P2 == ch1[tInd])[0]
                if ch1[tInd] == 0:
                    ch1[tInd] = P2[i]
                    break
                else:
                    continue
    for i in range(P2.shape[0]):
        if ch1[i] == 0:
            ch1[i] = P2[i]

    for i in range(int(cInd), int(cInd2)):
        if P1[i] not in ch2:
            tInd = i
            while True:
                tInd = np.where(P1 == ch2[tInd])[0]
                if ch2[tInd] == 0:
                    ch2[tInd] = P1[i]
                    break
                else:
                    continue
    for i in range(P1.shape[0]):
        if ch2[i] == 0:
            ch2[i] = P1[i]

    return ch1, ch2


def combPerORX(P1, P2):
    """
    Ordered crossover for permutated values.
    """

    ch1 = np.zeros_like(P1)
    ch2 = np.zeros_like(P1)

    rng = np.random.default_rng()
    cInds = np.array([])
    while len(np.unique(cInds)) < 2:
        cInds = rng.integers(1, P1.shape[0]+1, size=2)
    cInd = np.sort(cInds)[0]
    cInd2 = np.sort(cInds)[1]
    for i in range(cInd, cInd2):
        ch1[i], ch2[i] = P1[i], P2[i]

    P2 = np.hstack((P2[cInd2:], P2[:cInd2]))
    ch1 = np.hstack((ch1[cInd2:], ch1[:cInd2]))

    for i in range(P2.shape[0]):
        tInd = i
        while P2[i] not in ch1:
            if ch1[tInd] == 0:
                ch1[tInd] = P2[i]
                break
            elif tInd == P2.shape[0]-1:
                tInd = 0
                continue
            else:
                tInd += 1
                continue

    ch1 = np.hstack((ch1[ch1.shape[0]-cInd2:], ch1[:ch1.shape[0]-cInd2]))

    P1 = np.hstack((P1[cInd2:], P1[:cInd2]))
    ch2 = np.hstack((ch2[cInd2:], ch2[:cInd2]))

    for i in range(P1.shape[0]):
        tInd = i
        while P1[i] not in ch2:
            if ch2[tInd] == 0:
                ch2[tInd] = P1[i]
                break
            elif tInd == P1.shape[0]-1:
                tInd = 0
                continue
            else:
                tInd += 1
                continue

    ch2 = np.hstack((ch2[ch2.shape[0]-cInd2:], ch2[:ch2.shape[0]-cInd2]))

    return ch1, ch2


# -----------------------------------------
# -----------------------------------------
# Mutation
# -----------------------------------------
# -----------------------------------------


def mutBinBFL(chro, pm=0.2):
    """
    Bitflipping.
    """

    num = int(np.ceil(len(chro)*pm))
    rng = np.random.default_rng()
    inds = np.array([])
    while len(np.unique(inds)) < num:
        inds = rng.integers(0, len(chro), size=num)

    for i in inds:
        chro[i] = 1 - chro[i]

    return chro


# ----------------------------------------


def mutIntRan(chro, lb=None, ub=None, pm=0.2):
    """
    Random integer mutation.
    """

    if lb == None:
        lb = int(min(chro))
    if ub == None:
        ub = int(max(chro))

    num = int(np.ceil(len(chro)*pm))
    rng = np.random.default_rng()

    inds = np.array([])
    while len(np.unique(inds)) < num:
        inds = rng.integers(0, len(chro), size=num)

    vals = rng.integers(lb, ub+1, size=num)

    for i, v in zip(inds, vals):
        chro[i] = v

    return chro


def mutIntCre(chro, steps=[-1, 1], prob=[0.5, 0.5], pm=0.2):
    """
    Creep mutation for integer values.
    """

    CDF = [sum(prob[0:i]) for i in range(1, len(prob)+1)]
    num = int(np.ceil(len(chro)*pm))
    rng = np.random.default_rng()
    inds = np.array([])
    while len(np.unique(inds)) < num:
        inds = rng.integers(0, len(chro), size=num)

    for ind in inds:
        crus = rng.random()
        for i in range(len(CDF)):
            if crus > CDF[i]:
                continue
            else:
                chro[ind] += steps[i]
                break

    return chro


# ----------------------------------------


def mutConUni(chro, lb=None, ub=None, pm=0.1):
    """
    Uniform mutation for real values.
    """

    if lb == None:
        lb = int(min(chro))
    if ub == None:
        ub = int(max(chro))

    num = int(np.ceil(len(chro)*pm))
    rng = np.random.default_rng()

    inds = np.array([])
    while len(np.unique(inds)) < num:
        inds = rng.integers(0, len(chro), size=num)

    vals = (ub-lb)*rng.random(size=num)+lb

    for i, v in zip(inds, vals):
        chro[i] = v

    return chro


def mutConGau(chro, scale=1.0, pm=0.2):
    """
    Gaussian distribution mutation for real values.
    """

    num = int(np.ceil(len(chro)*pm))
    rng = np.random.default_rng()

    inds = np.array([])
    while len(np.unique(inds)) < num:
        inds = rng.integers(0, len(chro), size=num)

    for i in inds:
        chro[i] = rng.normal(chro[i], scale)

    return chro


# ----------------------------------------


def mutPerSwa(chro):
    """
    Swap mutation for permutated values.
    """

    rng = np.random.default_rng()
    inds = np.array([])
    while len(np.unique(inds)) < 2:
        inds = rng.integers(0, len(chro), size=2)

    chro[inds[0]], chro[inds[1]] = chro[inds[1]], chro[inds[0]]
    return chro


def mutPerIns(chro):
    """
    Insert mutation.
    """

    rng = np.random.default_rng()
    inds = np.array([])
    while len(np.unique(inds)) < 2:
        inds = rng.integers(0, len(chro)-1, size=2)
    inds = np.sort(inds)

    val = chro[inds[1]]
    chro = np.insert(np.delete(chro, inds[1]), inds[0]+1, val)

    return chro


def mutPerScr(chro):
    """
    Scramble mutation.
    """

    rng = np.random.default_rng()
    inds = np.array([])
    while len(np.unique(inds)) < 2:
        inds = rng.integers(0, len(chro), size=2)
    inds = np.sort(inds)

    chro[inds[0]:inds[1]+1] = rng.choice(chro[inds[0]:inds[1]+1],
                                         size=len(chro[inds[0]:inds[1]+1]), replace=False)

    return chro


def mutPerInv(chro):
    """
    Inverse mutation.
    """

    rng = np.random.default_rng()
    inds = np.array([])
    while len(np.unique(inds)) < 2:
        inds = rng.integers(0, len(chro), size=2)
    inds = np.sort(inds)

    chro[inds[0]:inds[1]+1] = chro[inds[0]:inds[1]+1][::-1]

    return chro


# -----------------------------------------
# -----------------------------------------
# Survivor Selection
# -----------------------------------------
# -----------------------------------------


def surSelAge(old, childs):
    """
    Age based survivor selection.
    """

    if old.shape[0] > childs.shape[0]:
        n = childs.shape[0]
        sur = np.delete(old, range(n), 0)
        sur = np.append(sur, childs, 0)
    else:
        sur = childs.copy()

    return sur


def surSelRep(old, childs, oldFit, chiFit, dup=False):
    """
    Replace worst survivor selection.
    """

    n = old.shape[0]

    fits = oldFit + chiFit
    fits = np.array(fitScale(fits))
    popu = np.append(old, childs, axis=0)

    if dup == False:
        popu, inds = np.unique(popu, axis=0, return_index=True)
        fits = fits[inds]

    args = np.argsort(0-fits)
    # mask = args < n
    sur = popu[args[0:n]]
    return sur
