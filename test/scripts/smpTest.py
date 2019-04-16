# python 3.7
# launch runEXE.py launch switchExp.exe
import subprocess
import sys
import os

#executable          =       'D:/GitHub/Repo/SMP/x64/Release/SMP_1271.exe'
#data                =       'D:/GitHub/Repo/SMP/test/data/part_graph_5.txt'
# formulation         =       '4'
# callbackOption      =       '3'
# for s in sys.argv:
#     if s.index == 0:
#         continue

dataLocation_1      =       sys.argv[1]  #ex: random_graph
dataLocation_2      =       sys.argv[2]  #ex: plan_graph
dataLocation_3      =       sys.argv[3]  #ex: group_1   
dataLocation_4      =       sys.argv[4]  #ex: dataset_1_1_1
samplesBit          =       sys.argv[5]  #ex: 1023(first ten graphs)
formulation         =       sys.argv[6]  #1\2\3\4
callbackOption      =       sys.argv[7]  #0\1\2\3
relaxOption         =       sys.argv[8]  #0\1
ns_sep_opt          =       sys.argv[9]  #0\1
timeLimit           =       sys.argv[10] #ex: 3600
cutNumber           =       sys.argv[11] #ex: 5
epsilon             =       sys.argv[12] #ex: 0.2       


os.chdir('../..') # to ...SMP/
#print(os.getcwd())
cwd = os.getcwd()
exeAbsltLocation    =       cwd + '\\x64\\Release\\SMP_1271.exe'
dataAbsltLocation   =       cwd + '\\test\\data\\'

#D:~/SMP/test/data/random_graph/plan_random/group_1/dataset1_1_1_2/
dataAbsltLocation = dataAbsltLocation + dataLocation_1 + '\\' + dataLocation_2 + '\\' + dataLocation_3 + '\\' + dataLocation_4 + '\\'

# if sys.argv[13] == '-s':
#     tempDataLocation = dataAbsltLocation + 'animal_' + sys.argv[14] + '.txt'
#     print ('\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
#               ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
#     #print ('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
#     print ('graph_'+sys.argv[14]+ ' START')
#     #print (tempDataLocation)
#     #subprocess.Popen([executable, tempDataLocation, formulation, callbackOption, relaxOption, ns_sep_opt, timeLimit, cutNumber, epsilon]).wait()
#     subprocess.Popen([exeAbsltLocation, tempDataLocation, formulation, callbackOption, relaxOption, ns_sep_opt, timeLimit, cutNumber, epsilon]).wait()
#     print ('graph_'+sys.argv[14]+' DONE')
    
# else:
#     for i in range(1, int(number_sample)+1):
#         tempDataLocation = ''
#         #D:/GitHub/Repo/SMPtest/data/random_graph/plan_random/group_1/dataset1_1_1_2/animal_10_2_5_84%_
#         tempDataLocation = dataAbsltLocation + 'animal_' + str(i) + '.txt'
#         print ('\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
#                 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
#         #print ('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
#         print ('graph_'+str(i)+ ' START')
#         #print (tempDataLocation)
#         #subprocess.Popen([executable, tempDataLocation, formulation, callbackOption, relaxOption, ns_sep_opt, timeLimit, cutNumber, epsilon]).wait()
#         subprocess.Popen([exeAbsltLocation, tempDataLocation, formulation, callbackOption, relaxOption, ns_sep_opt, timeLimit, cutNumber, epsilon]).wait()
#         print ('graph_'+str(i)+' DONE')

i = 1
for b in reversed(bin(int(samplesBit))):
    if b == 'b':
        i -= 2
        break
    elif b == str(1):
        tempDataLocation = ''
        #D:/GitHub/Repo/SMPtest/data/random_graph/plan_random/group_1/dataset1_1_1_2/animal_10_2_5_84%_
        tempDataLocation = dataAbsltLocation + 'animal_' + str(i) + '.txt'
        print ('\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
            ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
        #print ('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
        print ('graph_'+str(i)+ ' START')
        #print (tempDataLocation)
        #subprocess.Popen([executable, tempDataLocation, formulation, callbackOption, relaxOption, ns_sep_opt, timeLimit, cutNumber, epsilon]).wait()
        subprocess.Popen([exeAbsltLocation, tempDataLocation, formulation, callbackOption, relaxOption, ns_sep_opt, timeLimit, cutNumber, epsilon]).wait()
        print ('graph_'+str(i)+' DONE')
    i += 1    

#print('i = ', str(i))   

##############################
#### Begin to analyse data ###
##############################
from openpyxl import Workbook

if formulation == '1':
    result = open(dataAbsltLocation + "1_SCF.txt","r")
elif formulation == '2':
    result = open(dataAbsltLocation + "1_MCF.txt","r")
elif formulation == '3':
    result = open(dataAbsltLocation + "1_STEINER.txt","r")
elif formulation == '4':
    result = open(dataAbsltLocation + "1_NS.txt","r")

for line in result:
    print(line)

wb = Workbook()
ws = wb.active
ws.title = "animal_data"

