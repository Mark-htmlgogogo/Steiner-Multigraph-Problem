# python 3.7
# launch runEXE.py launch switchExp.exe
import subprocess
import sys

#executable          =       'D:/GitHub/Repo/SMP/x64/Release/SMP_1271.exe'
#data                =       'D:/GitHub/Repo/SMP/test/data/part_graph_5.txt'
# formulation         =       '4'
# callbackOption      =       '3'

executable          =       sys.argv[1]  #ex: SMP_1271.exe
dataLocation_1      =       sys.argv[2]  #ex: random_graph
dataLocation_2      =       sys.argv[3]  #ex: plan_graph
dataLocation_3      =       sys.argv[4]  #ex: group_1   
dataLocation_4      =       sys.argv[5]  #ex: dataset_1_1_1
dataFileName        =       sys.argv[6]  #ex: animal_10_2_5_84%_
number_sample       =       sys.argv[7]  #ex: 10
formulation         =       sys.argv[8]  #1\2\3\4
callbackOption      =       sys.argv[8]  #0\1\2\3
relaxOption         =       sys.argv[10] #0\1
ns_sep_opt          =       sys.argv[11] #0\1
timeLimit           =       sys.argv[12] #ex: 3600
cutNumber           =       sys.argv[13] #ex: 5
epsilon             =       sys.argv[14] #ex: 0.2

exeAbsltLocation    =       'D:\\aMain\\git\\x64\\Release\\'
dataAbsltLocation   =       'D:\\aMain\\git\\test\\data\\'

#D:/GitHub/Repo/SMPtest/data/random_graph/plan_random/group_1/dataset1_1_1_2/animal_10_2_5_84%_
exeAbsltLocation  = exeAbsltLocation + executable
dataAbsltLocation = dataAbsltLocation + dataLocation_1 + '\\' + dataLocation_2 + '\\' + dataLocation_3 + '\\' + dataLocation_4 + '\\' + dataFileName

for i in range(1, int(number_sample)+1):
    tempDataLocation = ''
    #D:/GitHub/Repo/SMPtest/data/random_graph/plan_random/group_1/dataset1_1_1_2/animal_10_2_5_84%_
    tempDataLocation = dataAbsltLocation + str(i) + '.txt'
    print ('\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
              ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
    #print ('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
    print ('graph_'+str(i)+ ' START')
    #print (tempDataLocation)
    #subprocess.Popen([executable, tempDataLocation, formulation, callbackOption, relaxOption, ns_sep_opt, timeLimit, cutNumber, epsilon]).wait()
    p = subprocess.Popen([exeAbsltLocation, tempDataLocation, formulation, callbackOption, relaxOption, ns_sep_opt, timeLimit, cutNumber, epsilon])
    p.wait()
    print ('graph_'+str(i)+' DONE')
    #print (tempDataLocation)

#p = subprocess.Popen([executable, data, formulation, callbackOption, relaxOption], stdout=subprocess.PIPE, 
#                        stderr=subprocess.PIPE, universal_newlines=True)
#process = subprocess.Popen([executable, data, formulation, callbackOption])
#out = p.communicate()
#print(err)
#print(out+"OUTOUTOUTOUTOUT")
#for line in out : print(line)
#print(out, end='\n')
#errcode = process.returncode
#subprocess.call([executable, data, formulation, callbackOption])
#subprocess.call(['D:/GitHub/Repo/SMP/x64/Release/SMP_1271.exe', 'D:/GitHub/Repo/SMP/test/data/part_graph_5.txt', '4', '3'])
#ubprocess.call([shutil.which("helloWorld.exe")])