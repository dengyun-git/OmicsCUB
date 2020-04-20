from scipy.stats import multinomial
import numpy


def Efor(cLeng, subLeng, cfPart)

    P = cfPart

    X = subLeng * cfPart

    Xtest = mod(X,1) * 10

    if Xtest.search[~0]=0:

        Pmax = multinomial.pmf(Xtest, sum(Xtest), P)

    else:

        Xpre = X

        rmv = subLeng - sum(np.floor(Xpre))

        Pmax = getPmax(cLeng, cfPart, np.floor(Xpre) , rmv)

    return Pmax
