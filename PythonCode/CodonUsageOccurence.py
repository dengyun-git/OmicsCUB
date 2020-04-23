import sys
import os
from scipy.stats import multinomial


##### function to get codon occurrence configruation

def getCOC(file1, file2, aaID, globalR):

    CodonTable = [['GAT', 'GAC'], ['ATA', 'ATT', 'ATC']]

    aaList = ['E', 'I']

    codonChoice = CodonTable[aaID]

    choiceSize = len(codonChoice)

    f1 = open(file1, 'r')

    while(True):

        line = f1.readline()

        if line == '':
            break

        if line[0] != '>':

            DNA = line
            Seq = []

            for i in range(len(DNA)-2):

                Seq.append(DNA[i*3: (i+1)*3])

            X = []

            for k in range(choiceSize):

                X.append(DNA.count(codonChoice[k]))

            lenSub = sum(X)

            Pi = multinomial.pmf(X, lenSub, globalR)

            Pmax = Efor(choiceSize, lenSub, globalR)

            f2 = open(file2, 'a')

            f2.write(aaList[aaID] + ',' + X + ',' + lenSub + ',' + Pi + ',' + Pmax \
                     + '\n')

    f2.close()
    f1.close()

    return
