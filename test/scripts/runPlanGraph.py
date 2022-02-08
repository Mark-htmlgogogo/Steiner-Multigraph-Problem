# Usage:
# py smpTest_smallscale.py  random_graph plan_random group_1 tg 19 4 3 0  0 3600 1200 0.2 1200 0.2
import subprocess
import sys
import os
from analyzePlanGraph import analyze, analyzeLB, writeStartIndex

dataLocation_1 = sys.argv[1]  # ex: random_graph
dataLocation_2 = sys.argv[2]  # ex: plan_graph
dataLocation_3 = sys.argv[3]  # ex: group_1
dataLocation_4 = sys.argv[4]  # ex: dataset_1_1_1
graph_number = sys.argv[5]  # ex: 1023(first ten graphs)
formulation = sys.argv[6]  # 1\2\3\4
callback_option = sys.argv[7]  # 0\1\2\3
relax_option = sys.argv[8]  # 0\1
ns_sep_opt = sys.argv[9]  # 0\1
LB_MaxRestart = sys.argv[10]
LB_MaxIter = sys.argv[11]
Rmin = sys.argv[12]
Rmax = sys.argv[13]
BCSolNum = sys.argv[14]
BCTime = sys.argv[15]
MIPDisplayLevel = sys.argv[16]
time_limit = sys.argv[17]  # ex: 3600
max_cut_number_lazy = sys.argv[18]  # ex: 5
epsilon_lazy = sys.argv[19]  # ex: 0.2
max_cut_number_user = sys.argv[20]  # ex: 5
epsilon_user = sys.argv[21]  # ex: 0.2
UseLocalBranch = sys.argv[22]
LB_CP_Option = sys.argv[23]
lazy_sep_opt = sys.argv[24]

os.chdir('../..')  # to ...SMP/
cwd = os.getcwd()
exeAbsltLocation = cwd + '\\x64\\Release\\SMP_1271_test_ns.exe'
dataAbsltLocation = cwd + '\\test\\data\\'

# D:~/SMP/test/data/random_graph/plan_random/group_1/dataset1_1_1_2/
ndataAbsltLocation = dataAbsltLocation + dataLocation_1 + '\\' + \
    dataLocation_2 + '\\' + dataLocation_3 + '\\'

prefix = 'dataset1_1_'
#suffixPool = ['t_g']
suffixPool = ['5_1', '5_2', '5_3', '5_4']
# suffixPool = ['1_1', '1_2', '1_3', '1_4',
#               '2_1', '2_2', '2_3', '2_4',
#               '3_1', '3_2', '3_3', '3_4',
#               '4_1', '4_2', '4_3', '4_4',
#               '5_1', '5_2', '5_3', '5_4']
for suffix in suffixPool:
    fileLocation = ndataAbsltLocation + prefix + suffix + '\\'
    writeStartIndex(fileLocation, formulation)
    for i in range(1, int(graph_number)+1):
        tempDataLocation = fileLocation + 'animal_' + str(i) + '.txt'
        print('\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
            ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
        print(fileLocation+'\ngraph_'+str(i) + ' START')
        subprocess.Popen([exeAbsltLocation, tempDataLocation, formulation, callback_option, relax_option,
                          ns_sep_opt, LB_MaxRestart, LB_MaxIter, Rmin, Rmax, BCSolNum, BCTime, MIPDisplayLevel,
                          time_limit, max_cut_number_lazy, epsilon_lazy, max_cut_number_user, epsilon_user, UseLocalBranch,
                          LB_CP_Option, lazy_sep_opt]).wait()
        print(fileLocation+'\ngraph_'+str(i) + ' DONE')
    if int(UseLocalBranch) == 0:
        analyze(graph_number, formulation, fileLocation)
    else:
        analyzeLB(graph_number, formulation, fileLocation)
pass
