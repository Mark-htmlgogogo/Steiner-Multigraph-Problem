import subprocess
import sys
import os

dataLocation_1 = sys.argv[1]  # ex: random_graph
dataLocation_2 = sys.argv[2]  # ex: plan_graph
dataLocation_3 = sys.argv[3]  # ex: group_1
dataLocation_4 = sys.argv[4]  # ex: dataset_1_1_1
os.chdir('../..')  # to ...SMP/
cwd = os.getcwd()
exeAbsltLocation = cwd + '\\x64\\Release\\SMP_1271_test_ns.exe'
dataAbsltLocation = cwd + '\\test\\data\\'

# D:~/SMP/test/data/random_graph/plan_random/group_1/dataset1_1_1_2/
dataAbsltLocation1 = dataAbsltLocation + dataLocation_1 + '\\' + \
    dataLocation_2 + '\\' + dataLocation_3 + '\\' + \
    dataLocation_4 + '\\' + "1_NS.txt"

dataAbsltLocation2 = dataAbsltLocation + dataLocation_1 + '\\' + \
    dataLocation_2 + '\\' + dataLocation_3 + '\\' + \
    dataLocation_4 + '\\' + "1UserIndex.txt"

f1 = open(str(dataAbsltLocation1), "r")

line1 = f1.readlines()
f1.close()

epsilon_pool = ["0.001", "0.1", "0.2", "0.4", "0.5"]
cutnumber = ["1", "3", "5", "7", "-1"]
epsilon_lazy = "0.0"
max_cut_number_lazy = "-1"

strlazy = "max_cut_number_lazy: -1, epsilon_lazy: 0.0 \n "

idx = 0
graph_number = 50

f2 = open(str(dataAbsltLocation2), "w")
f2.write(strlazy + "\n")

for max_cut_number_user in cutnumber:
    for epsilon_user in epsilon_pool:

        usernumber = "max_cut_number_user: " + max_cut_number_user
        userepsilon = ", epailon_user: " + epsilon_user
        strpre = usernumber + userepsilon + "\n"

        AveTime = 0.0
        Nodes = 0
        Cuts = 0

        for i in range(1, graph_number+1):
            str_list = line1[50 * idx + i - 1].split()
            AveTime += (float(str_list[0]) / float(graph_number))
            Nodes += (float(str_list[1]) / float(graph_number))
            Cuts += (float(str_list[2]) / float(graph_number))

        _str1 = "Ave Time: " + str(round(AveTime, 3))
        _str1 += ",  Nodes: " + str(int(Nodes))
        _str1 += ",  Cuts:  " + str(int(Cuts)) + "\n"

        fstr = strpre + _str1
        print(fstr)
        f2.write(fstr)
        fstr = str(round(AveTime, 3))+"/" + \
            str(int(Nodes)) + "/" + str(int(Cuts)) + "\n\n"
        f2.write(fstr)
        idx += 1

f2.close()
print("----------------------Finsih-----------------------")
