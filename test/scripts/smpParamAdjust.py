# python 3.7
# launch runEXE.py launch switchExp.exe
from hyperopt import fmin, tpe, hp
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
dataAbsltLocation = dataAbsltLocation + dataLocation_1 + '\\' + \
    dataLocation_2 + '\\' + dataLocation_3 + '\\' + dataLocation_4 + '\\'


space = {
    'ns_sep_opt': hp.randint('ns_sep_opt', 1),
    'max_cut_number_lazy': hp.randint('max_cut_number_lazy', 1200)+1,
    'epsilon_lazy': hp.quniform('epsilon_lazy', 0.001, 1, 0.01),
    'max_cut_number_user': hp.randint('max_cut_number_user', 1200)+1,
    'epsilon_user': hp.quniform('epsilon_user', 0.001, 1, 0.01)
}


read_position = 0


def smp(space):
    ns_sep_opt = str(space['ns_sep_opt'])
    max_cut_number_lazy = str(space['max_cut_number_lazy'])
    epsilon_lazy = str(space['epsilon_lazy'])
    max_cut_number_user = str(space['max_cut_number_user'])
    epsilon_user = str(space['epsilon_user'])

    #########################
    ### run .exe in loop ####
    #########################
    number_trials = 1
    for b in reversed(bin(int(sampples_bit))):
        if b == 'b':
            number_trials -= 1
            break
        elif b == str(1):
            tempDataLocation = ''
            # D:/GitHub/Repo/SMPtest/data/random_graph/plan_random/group_1/dataset1_1_1_2/animal_10_2_5_84%_
            tempDataLocation = dataAbsltLocation + \
                'animal_' + str(number_trials) + '.txt'
            print('\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
                ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
            print('graph_'+str(number_trials) + ' START')
            subprocess.Popen([exeAbsltLocation, tempDataLocation, formulation, callback_option, relax_option,
                              ns_sep_opt, LB_MaxRestart, LB_MaxIter, Rmin, Rmax, BCSolNum, BCTime, MIPDisplayLevel,
                              time_limit, max_cut_number_lazy, epsilon_lazy, max_cut_number_user, epsilon_user,
                              UseLocalBranch, LB_CP_Option, lazy_sep_opt]).wait()
            print('graph_'+str(number_trials)+' DONE')
        number_trials += 1

    ##############################
    ### Begin to analyse data ####
    ##############################
    if formulation == '1':
        result_file_name = dataAbsltLocation + "1_SCF.txt"
        result_file_stream = open(result_file_name, "r")
    elif formulation == '2':
        result_file_name = dataAbsltLocation + "1_MCF.txt"
        result_file_stream = open(result_file_name, "r")
    elif formulation == '3':
        result_file_name = dataAbsltLocation + "1_STEINER.txt"
        result_file_stream = open(result_file_name, "r")
    elif formulation == '4':
        result_file_name = dataAbsltLocation + "1_NS.txt"
        result_file_stream = open(result_file_name, "r")

    total_time = 0
    argv_time = 0
    global read_position
    result_file_stream.seek(read_position, 0)
    for line in result_file_stream:
        split_result_line = line.split()
        total_time += float(split_result_line[1])
    argv_time = total_time / number_trials
    print('argv_time = ', argv_time)
    read_position = result_file_stream.tell()
    # open(result_file_name,'w').close() # empty the file.

    return argv_time


best = fmin(
    fn=smp,
    space=space,
    algo=tpe.suggest,
    max_evals=int(number_of_evals)
)

print(best)
