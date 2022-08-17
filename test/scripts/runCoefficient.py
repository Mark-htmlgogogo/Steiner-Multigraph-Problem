# Usage:
# run NS test with general random graph
import subprocess
import sys
import os
from datetime import datetime
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

subfolder = ["n1000_t30_p3_b0005_v06",
             "n2000_t30_p3_b0005_v06",
             "n2500_t30_p3_b0005_v06",
             "n3000_t30_p3_b0005_v06"]

coefficients = ["0.0", "0.2", "0.4", "0.6", "0.8", "1.0"]
# coefficients = ["0.0"]

dataAbsltLocation = cwd + '\\test\\data\\'
dataAbsltLocation = dataAbsltLocation + dataLocation_1 + '\\' + \
    dataLocation_2 + '\\' + dataLocation_3 + '\\' + dataLocation_4 + '\\'

# tmp = dataAbsltLocation
# for foder in subfolder:
#     dataAbsltLocation = tmp + foder + '\\'
myList = os.listdir(dataAbsltLocation)
myList.sort()

for file in myList:
    if file == "1_NS.txt" or file == "1_NS_relax.txt" or file == "1_MCF.txt":
        continue

    for coeff in coefficients:
        if file < '4-cc12-2nFRS.txt':
            continue
        elif file > '4-cc12-2nFRS.txt':
            pass
        elif coeff < '0.1':
            continue

        if relax_option == "0":
            f = open(dataAbsltLocation+"1_NS.txt", "a")
        elif relax_option == "1":
            f = open(dataAbsltLocation+"1_NS_relax.txt", "a")

        f.write("\n\n"+str(coeff)+"\n")
        f.flush()
        f.close()

        tempDataLocation = dataAbsltLocation + file
        print('\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
                    ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
        print(file + ' START' + " at ")

        now = datetime.now()
        current_time = now.strftime("%H:%M:%S")
        print("Current Time =", current_time)

        subprocess.Popen([exeAbsltLocation, tempDataLocation, formulation, callback_option, relax_option,
                          ns_sep_opt, LB_MaxRestart, LB_MaxIter, Rmin, Rmax, BCSolNum, BCTime, MIPDisplayLevel,
                          time_limit, max_cut_number_lazy, epsilon_lazy, max_cut_number_user, epsilon_user, UseLocalBranch, LB_CP_Option, lazy_sep_opt, coeff]).wait()
        print(file + ' DONE')
   # get_last_line(len(list(filter(lambda x: x[0] == 'a', myList))),formulation, UseLocalBranch, dataAbsltLocation)
