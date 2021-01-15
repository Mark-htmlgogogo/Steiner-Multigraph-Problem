import os
import time
import datetime
import math
import string

formulationResultPool = ["1_SCF.txt",
                         "1_MCF.txt", "1_STEINER.txt", "1_NS.txt"]


def writeStartIndex(dataAbsltLocation, formulation):
    dataStoreFile = dataAbsltLocation + \
        formulationResultPool[int(formulation)-1]
    print(dataStoreFile)
    f1 = open(str(dataStoreFile), "a")
    f1.write('#\n')
    f1.close()
    pass


def analyze(graph_number, formulation, dataAbsltLocation):
    dataStoreFile = dataAbsltLocation + \
        formulationResultPool[int(formulation)-1]
    f1 = open(str(dataStoreFile), "r")
    lines = f1.readlines()
    newLine = []
    for i in range(1, int(graph_number)+1):
        str_list = lines[-i].split()
        if len(str_list) == 0:
            continue
        elif len(str_list) == 1 and str_list[0] == '#':
            break
        else:
            newLine.append(str_list)

    # * start compose data
    D = [[0.0]*8 for i in range(2)]
    D[0][6] = D[1][6] = int(0)
    for line in newLine:
        # * param define
        time = float(line[1])
        finalGap = float(line[4])
        pos = int(0)

        if finalGap > 1e-5:
            pos = 1

        D[pos][0] += 1
        D[pos][1] += time
        D[pos][2] += finalGap
        D[pos][3] = max(D[pos][3], finalGap)

    _size = len(newLine)
    diffGraphNumber = int(graph_number) - _size
    pos = 1

    D[pos][0] += diffGraphNumber
    D[pos][1] += diffGraphNumber*3600.0
    D[pos][2] += diffGraphNumber*1
    D[pos][3] = 1.0
    f1.close()

    # * write to file
    f1 = open(str(dataStoreFile), "a")
    s = ''

    s += str(D[0][0])
    for i in range(3):
        if i == 0:
            continue
        s += '   '+str(round(D[0][i]/float(D[0][0]), 5))
    s += '   '+str(round(D[0][3], 5))

    if D[1][0] != 0:
        s += '\n' + str(D[1][0])
        for i in range(3):
            if i == 0:
                continue
            s += '   '+str(round(D[1][i]/float(D[1][0]), 5))
        s += '   '+str(round(D[1][3], 5))

    f1.write('================Analyze Result===============\n')
    f1.write(s)
    f1.write('\n=============================================')
    pass


def analyzeLB():
    pass


if __name__ == "__main__":
    testLocation = 'D:\\aMain\\git\\test\\scripts\\'
    #writeStartIndex(testLocation, '4')
    analyze(5, 4, testLocation)
    pass
