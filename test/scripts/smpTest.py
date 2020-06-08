import subprocess
import sys
import os

dataLocation_1 = sys.argv[1]  # ex: random_graph
dataLocation_2 = sys.argv[2]  # ex: plan_graph
dataLocation_3 = sys.argv[3]  # ex: group_1
dataLocation_4 = sys.argv[4]  # ex: dataset_1_1_1
sampples_bit = sys.argv[5]  # ex: 1023(first ten graphs)
formulation = sys.argv[6]  # 1\2\3\4
callback_option = sys.argv[7]  # 0\1\2\3
relax_option = sys.argv[8]  # 0\1
ns_sep_opt = sys.argv[9]  # 0\1
time_limit = sys.argv[10]  # ex: 3600
max_cut_number_lazy = sys.argv[11]  # ex: 5
epsilon_lazy = sys.argv[12]  # ex: 0.2
max_cut_number_user = sys.argv[13]  # ex: 5
epsilon_user = sys.argv[14]  # ex: 0.2

os.chdir('../..')  # to ...SMP/
cwd = os.getcwd()
exeAbsltLocation = cwd + '\\x64\\Release\\SMP_1271_test_ns.exe'
dataAbsltLocation = cwd + '\\test\\data\\'

# D:~/SMP/test/data/random_graph/plan_random/group_1/dataset1_1_1_2/
dataAbsltLocation = dataAbsltLocation + dataLocation_1 + '\\' + \
    dataLocation_2 + '\\' + dataLocation_3 + '\\' + dataLocation_4 + '\\'

i = 1
for b in reversed(bin(int(sampples_bit))):
    if b == 'b':
        i -= 2
        break
    elif b == str(1):
        tempDataLocation = ''
        # D:/GitHub/Repo/SMPtest/data/random_graph/plan_random/group_1/dataset1_1_1_2/animal_10_2_5_84%_
        tempDataLocation = dataAbsltLocation + 'animal_' + str(i) + '.txt'
        print('\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
            ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
        print('graph_'+str(i) + ' START')
        subprocess.Popen([exeAbsltLocation, tempDataLocation, formulation, callback_option, relax_option,
                          ns_sep_opt, time_limit, max_cut_number_lazy, epsilon_lazy, max_cut_number_user, epsilon_user]).wait()
        print('graph_'+str(i)+' DONE')
    i += 1

