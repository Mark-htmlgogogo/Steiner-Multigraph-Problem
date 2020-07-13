# Usage:
# py smpTest_smallscale.py  random_graph plan_random group_1 tg 19 4 3 0  0 3600 1200 0.2 1200 0.2
import subprocess
import sys
import os
from data_calculate import get_last_line

dataLocation_1 = sys.argv[1]  # ex: random_graph
dataLocation_2 = sys.argv[2]  # ex: plan_graph
dataLocation_3 = sys.argv[3]  # ex: group_1
dataLocation_4 = sys.argv[4]  # ex: dataset_1_1_1
st_file = sys.argv[5]
ed_file = sys.argv[6]
graph_number = sys.argv[7]  # ex: 1023(first ten graphs)
formulation = sys.argv[8]  # 1\2\3\4
callback_option = sys.argv[9]  # 0\1\2\3
relax_option = sys.argv[10]  # 0\1
ns_sep_opt = sys.argv[11]  # 0\1
LB_MaxRestart = sys.argv[12]
LB_MaxIter = sys.argv[13]
Rmin = sys.argv[14]
Rmax = sys.argv[15]
BCSolNum = sys.argv[16]
BCTime = sys.argv[17]
MIPDisplayLevel = sys.argv[18]
time_limit = sys.argv[19]  # ex: 3600
max_cut_number_lazy = sys.argv[20]  # ex: 5
epsilon_lazy = sys.argv[21]  # ex: 0.2
max_cut_number_user = sys.argv[22]  # ex: 5
epsilon_user = sys.argv[23]  # ex: 0.2
UseLocalBranch = sys.argv[24]
LB_CP_Option = sys.argv[25]
lazy_sep_opt = sys.argv[26]

os.chdir('../..')  # to ...SMP/
cwd = os.getcwd()
exeAbsltLocation = cwd + '\\x64\\Release\\SMP_1271_test_ns.exe'
dataAbsltLocation = cwd + '\\test\\data\\'

runformulation = ["4", "2"]

for idx in range(int(st_file), int(ed_file) + 1):
    tmp = dataLocation_4[0:len(dataLocation_4) - 1]
    dataLocation_4 = tmp + str(idx)

    for formulation_type in runformulation:
        formulation = str(formulation_type)

        dataAbsltLocation = cwd + '\\test\\data\\'
        dataAbsltLocation = dataAbsltLocation + dataLocation_1 + '\\' + \
            dataLocation_2 + '\\' + dataLocation_3 + '\\' + dataLocation_4 + '\\'

        for i in range(2, int(graph_number)+1):
            tempDataLocation = ''
            # D:/GitHub/Repo/SMPtest/data/random_graph/plan_random/group_1/dataset1_1_1_2/animal_10_2_5_84%_
            tempDataLocation = dataAbsltLocation + 'animal_' + str(i) + '.txt'
            print('\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
                    ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
            print('graph_'+str(i) + ' START')
            subprocess.Popen([exeAbsltLocation, tempDataLocation, formulation, callback_option, relax_option,
                              ns_sep_opt, LB_MaxRestart, LB_MaxIter, Rmin, Rmax, BCSolNum, BCTime, MIPDisplayLevel,
                              time_limit, max_cut_number_lazy, epsilon_lazy, max_cut_number_user, epsilon_user, UseLocalBranch,
                              LB_CP_Option, lazy_sep_opt]).wait()
            print('graph_'+str(i)+' DONE')

        # calculate data
        get_last_line(graph_number, formulation,
                      UseLocalBranch, dataAbsltLocation)
