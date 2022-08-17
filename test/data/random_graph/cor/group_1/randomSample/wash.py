import os
from tkinter import font
from numpy import random
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import skewnorm


dirs = "D:/aMain/git/test/data/random_graph/cor/group_1/randomSample"
files = []

# * 50 - 1000
x = np.linspace(50, 1000, 500)

# Varying positional arguments
a = 4
loc = 200
scale = 200

y1 = skewnorm .pdf(x, a, loc, scale)
plt.figure()
plt.rcParams.update({'font.size': 14})
plt.plot(x, y1, "*")
plt.xticks(np.linspace(50, 1000, 5))
plt.savefig(
    'D:/aMain/git/test/data/random_graph/cor/group_1/randomSample/skenorm.pdf', dpi=1200)

for file in files:
    filePath = dirs+'/'+file
    newFilePath = dirs+'/new_' + file
    fr = open(filePath, 'r')
    fw = open(newFilePath, 'w')

    l = fr.readline()
    fw .writelines(l)
    N = int(l.split()[0])
    for i in range(0, N):
        l = fr.readline().split()
        l[1] = str(int(skewnorm.rvs(a, loc, scale)))
        nl = l[0] + ' ' + l[1] + '\n'
        fw.writelines(nl)
    while True:
        l = fr.readline()
        fw.writelines(l)
        if not l:
            break
    fr.close()
    fw.close()
    print('finish'+file)
