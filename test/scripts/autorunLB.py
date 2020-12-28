# Usage:
# run LB within general random graph
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

prefixpool = ["3000"]
MaxIterPool = ["2","10", "20"]
LBTimeLimitPool = ["5","10","20","50"]
LB_MaxRestart = "1"
Rmin = "10"
Rmax = "30"
BCSolNum = "10"

# pre-define of the function



def get_last_line(graph_number, formulation, UseLocalBranch, dataAbsltLocation, LB_MaxRestart, LB_MaxIter, BCTime):
    # calculate data

    if formulation == "1":
        dataStoreFile = dataAbsltLocation + "1_SCF.txt"
    elif formulation == "2":
        dataStoreFile = dataAbsltLocation + "1_MCF.txt"
    elif formulation == "3":
        dataStoreFile = dataAbsltLocation + "1_STEINER.txt"
    elif formulation == "4":
        dataStoreFile = dataAbsltLocation + "1_NS.txt"

    Ave_TOT_TIME = 0.0
    Ave_LB_TIME = 0.0
    Ave_FINAL_SOLVE_TIME = 0.0
    Ave_LB_RUN_TIME = 0
    Final_solve_gap = 0.0
    Gap = 0.0
    Nodes = 0
    Cuts = 0

    f1 = open(str(dataStoreFile), "r")
    lines = f1.readlines()
    for i in range(1, int(graph_number)+1):
        str_list = lines[-i].split()
        if len(str_list) == 0:
            Ave_TOT_TIME += (float(3600)*(float(graph_number)-i+1) /
                             float(graph_number))
            #Ave_LB_TIME += (float(3600) / float(graph_number))
            Ave_FINAL_SOLVE_TIME += (float(3600) *
                                     (float(graph_number)-i+1) / float(graph_number))
            #Ave_LB_RUN_TIME += (float(str_list[4]) / float(graph_number))
            Final_solve_gap += (float(1)*(float(graph_number) -
                                          i+1) / float(graph_number))
            Gap += (float(1)*(float(graph_number)-i+1) / float(graph_number))
            #Nodes += (int(str_list[7]) / float(graph_number))
            #Cuts += (int(str_list[8]) / float(graph_number))
            break
        else:
            Ave_TOT_TIME += (float(str_list[1]) / float(graph_number))
            Ave_LB_TIME += (float(str_list[2]) / float(graph_number))
            Ave_FINAL_SOLVE_TIME += (float(str_list[3]) / float(graph_number))
            Ave_LB_RUN_TIME += (float(str_list[4]) / float(graph_number))
            Final_solve_gap += (float(str_list[5]) / float(graph_number))
            Gap += (float(str_list[6]) / float(graph_number))
            Nodes += (int(str_list[7]) / float(graph_number))
            Cuts += (int(str_list[8]) / float(graph_number))

    f1.close()

    f1 = open(str(dataStoreFile), "a")
    _str1 = "Ave_TOT_TIME = " + str(round(Ave_TOT_TIME, 5))
    _str1 += "\nAve_LB_TIME = " + str(round(Ave_LB_TIME, 5))
    _str1 += "\nAve_FINAL_SOLVE_TIME = " + \
        str(round(Ave_FINAL_SOLVE_TIME, 5))
    _str1 += "\nAve_LB_RUN_TIME = " + str(round(Ave_LB_RUN_TIME, 5))
    _str1 += "\nFinal solve gap = " + str(round(Final_solve_gap, 5))
    _str1 += "\nGaps = " + str(round(Gap, 5))
    _str1 += "\nNodes = " + str(int(Nodes))
    _str1 += "\nCuts = " + str(int(Cuts))
    _str1 += "\nMaxRestart = " + LB_MaxRestart + \
        "   MaxIter = " + LB_MaxIter + "   BCTime: " + BCTime
    print(_str1)
    f1.write("-----------------HERE----------------------\n")
    f1.write(_str1)
    f1.write("\n-----------------END----------------------\n\n")
    f1.close()


# end

# begin normal process
for idx in prefixpool:
    dataLocation_4 = 'n' + idx + "_t30_p3_b0005_v06"

    dataAbsltLocation = cwd + '\\test\\data\\'
    dataAbsltLocation = dataAbsltLocation + dataLocation_1 + '\\' + \
        dataLocation_2 + '\\' + dataLocation_3 + '\\' + dataLocation_4 + '\\'

    for BCTime in LBTimeLimitPool:
        for LB_MaxIter in MaxIterPool:
            for i in range(1, int(graph_number)+1):
                tempDataLocation = ''
                # D:/GitHub/Repo/SMPtest/data/random_graph/plan_random/group_1/dataset1_1_1_2/animal_10_2_5_84%_
                tempDataLocation = dataAbsltLocation + \
                    'animal_' + str(i) + '.txt'
                print('\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
                        ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
                print('graph_' + str(i) + ' START')
                subprocess.Popen([exeAbsltLocation, tempDataLocation, formulation, callback_option, relax_option,
                                  ns_sep_opt, LB_MaxRestart, LB_MaxIter, Rmin, Rmax, BCSolNum, BCTime, MIPDisplayLevel,
                                  time_limit, max_cut_number_lazy, epsilon_lazy, max_cut_number_user, epsilon_user, UseLocalBranch,
                                  LB_CP_Option, lazy_sep_opt]).wait()
                print('graph_' + str(i)+' DONE')

            # calculate data
            get_last_line(graph_number, formulation,
                          UseLocalBranch, dataAbsltLocation, LB_MaxRestart, LB_MaxIter, BCTime)
