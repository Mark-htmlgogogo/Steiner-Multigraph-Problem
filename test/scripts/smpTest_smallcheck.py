import subprocess
import sys
import os

dataLocation_1 = sys.argv[1]  # ex: random_graph
dataLocation_2 = sys.argv[2]  # ex: plan_graph
dataLocation_3 = sys.argv[3]  # ex: group_1
dataLocation_4 = sys.argv[4]  # ex: dataset_1_1_1
compare_file1 = sys.argv[5]  # 1_MCF
compare_file2 = sys.argv[6]  # 1_NS

os.chdir('../..')  # to ...SMP/
cwd = os.getcwd()
exeAbsltLocation = cwd + '\\x64\\Release\\SMP_1271_test_ns.exe'
dataAbsltLocation = cwd + '\\test\\data\\'

# D:~/SMP/test/data/random_graph/plan_random/group_1/dataset1_1_1_2/
dataAbsltLocation1 = dataAbsltLocation + dataLocation_1 + '\\' + \
    dataLocation_2 + '\\' + dataLocation_3 + '\\' + \
    dataLocation_4 + '\\' + compare_file1 + ".txt"

dataAbsltLocation2 = dataAbsltLocation + dataLocation_1 + '\\' + \
    dataLocation_2 + '\\' + dataLocation_3 + '\\' + \
    dataLocation_4 + '\\' + compare_file2 + ".txt"

f1 = open(str(dataAbsltLocation1), "r")
f2 = open(str(dataAbsltLocation2), "r")

line1 = f1.readlines()
f1.close()
line2 = f2.readlines()
f2.close()

count = 0
flag = 0
for i in range(len(line1)):
    if (int(line2[i]) != int(line1[i])):
        if (flag == 0):
            flag = 1
            print("Find dismatch:")
        print(i+1)
        count += 1

if flag == 0:
    print("No dismatch")
else:
    print("Dismatch number is", count)

print("----------------------Finsih-----------------------")
