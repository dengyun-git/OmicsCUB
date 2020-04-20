import os
import shutil

f = open('nameListABC.csv')

while(True):

    line = f.readline().strip('\n')

    if line == []:
        break

    else:
        fileName = line + '.fa'
        shutil.copy('/Users/yd46/Documents/MATLAB/' + fileName, '/Users/yd46/PycharmProjects/venv/')

f.close()
