# Usage:
# py smpTest_smallscale.py  random_graph plan_random group_1 tg 19 4 3 0  0 3600 1200 0.2 1200 0.2
import subprocess
import sys
import os

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


def get_last_line(graph_number, formulation, ns_sep_opt, dataAbsltLocation):
    # calculate data

    if formulation == "1":
        dataStoreFile = dataAbsltLocation + "1_SCF.txt"
    elif formulation == "2":
        dataStoreFile = dataAbsltLocation + "1_MCF.txt"
    elif formulation == "3":
        dataStoreFile = dataAbsltLocation + "1_STEINER.txt"
    elif formulation == "4":
        dataStoreFile = dataAbsltLocation + "1_NS.txt"

    AveTime = 0.0
    Nodes = 0
    Cuts = 0
    Gap = 0.0
    f1 = open(str(dataStoreFile), "r")
    lines = f1.readlines()
    for i in range(1, int(graph_number) + 1):
        # begin collect data
        str_list = lines[-i].split()
        AveTime += (float(str_list[1]) / float(graph_number))
        Nodes += (int(str_list[2]) / float(graph_number))
        Cuts += (int(str_list[3]) / float(graph_number))
        Gap += (float(str_list[4]) / float(graph_number))
    f1.close()
    f1 = open(str(dataStoreFile), "a")
    _str = " Ave Time: " + str(round(AveTime, 3))
    _str = _str + "  Nodes: " + str(int(Nodes))
    _str = _str + "  Cuts: " + str(int(Cuts))
    _str = _str + "  Gaps: " + str(round(Gap, 3))
    print(_str)
    f1.write("-----------------HERE----------------------\n")
    f1.write(_str)
    f1.write("\n-----------------END----------------------\n\n")
    f1.close()

    return


# D:~/SMP/test/data/random_graph/plan_random/group_1/dataset1_1_1_2/
dataAbsltLocation = dataAbsltLocation + dataLocation_1 + '\\' + \
    dataLocation_2 + '\\' + dataLocation_3 + '\\' + dataLocation_4 + '\\'
myList = os.listdir(dataAbsltLocation)
# myList.sort()
newList = []
for file in myList:
    if not file[0].isdigit():
        newList.append(file)
callBackPool = ["1", "3"]

for callback_option in callBackPool:
    if callback_option == "1":
        nsSepPool = ["1"]
    else:
        nsSepPool = ["0", "1"]
    for ns_sep_opt in nsSepPool:
        for file in newList:
            tempDataLocation = ''
            tempDataLocation = dataAbsltLocation + file
            print('\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
                        ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
            print(file + ' START')
            subprocess.Popen([exeAbsltLocation, tempDataLocation, formulation, callback_option, relax_option,
                              ns_sep_opt, LB_MaxRestart, LB_MaxIter, Rmin, Rmax, BCSolNum, BCTime, MIPDisplayLevel,
                              time_limit, max_cut_number_lazy, epsilon_lazy, max_cut_number_user, epsilon_user, UseLocalBranch,
                              LB_CP_Option, lazy_sep_opt]).wait()
            print(file + ' DONE')
        get_last_line(len(newList), formulation, ns_sep_opt, dataAbsltLocation)
