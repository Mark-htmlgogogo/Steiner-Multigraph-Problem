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

# 4 1 0 0 3 10 10 30 10 20 1 3600 1200 0.2 1200 0.2 1 0 1
# subfile = ["Copenhagen14", "MulPartition", "SPG-PUCN", "SPG-GAPS\\skutella", "SteinerLiB\\I320",
#            "SteinerLiB\\I640", "SteinerLiB\\MC", "SteinerLiB\\SP", "SteinerLiB\\X"]
# subfile = ["n1000_t30_p3_b0005_v06","n1500_t30_p3_b0005_v06","n2000_t30_p3_b0005_v06","n2500_t30_p3_b0005_v06","n3000_t30_p3_b0005_v06",
#            "n4000_t30_p3_b0005_v06","n5000_t30_p3_b0005_v06"]
#subfile = ["MulPartition"]
subfile = ["SPG-PUCN", "SPG-GAPS\\skutella", "SteinerLiB\\I320",
           "SteinerLiB\\I640", "SteinerLiB\\MC", "SteinerLiB\\SP", "SteinerLiB\\X"]

for doc in subfile:
    dataAbsltLocation = ''
    dataAbsltLocation = ndataAbsltLocation + doc + '\\'
    #dataAbsltLocation = ndataAbsltLocation
    myList = os.listdir(dataAbsltLocation)
    for file in myList:
        if file == "1_MCF.txt" or file == "1_NS.txt" or file == "0GraphInfo.txt" or file == "1_STEINER.txt":
            continue
        if doc == "SPG-PUCN" and file <= "cc3-11nFRS.txt":
            continue
        tempDataLocation = ''
        # D:/GitHub/Repo/SMPtest/data/random_graph/plan_random/group_1/dataset1_1_1_2/animal_10_2_5_84%_
        tempDataLocation = dataAbsltLocation + file
        print('\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
                ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
        print(file + ' START')
        subprocess.Popen([exeAbsltLocation, tempDataLocation, formulation, callback_option, relax_option,
                          ns_sep_opt, LB_MaxRestart, LB_MaxIter, Rmin, Rmax, BCSolNum, BCTime, MIPDisplayLevel,
                          time_limit, max_cut_number_lazy, epsilon_lazy, max_cut_number_user, epsilon_user, UseLocalBranch,
                          LB_CP_Option, lazy_sep_opt]).wait()
        print(file + ' DONE')
    # calculate data
    # get_last_line(str(len(myList)-3), formulation, UseLocalBranch, dataAbsltLocation)
