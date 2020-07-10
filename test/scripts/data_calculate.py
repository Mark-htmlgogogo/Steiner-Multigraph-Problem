import os
import time
import datetime
import math
import string


def get_last_line(graph_number, formulation, UseLocalBranch, dataAbsltLocation):
    # calculate data

    if formulation == "1":
        dataStoreFile = dataAbsltLocation + "1_SCF.txt"
    elif formulation == "2":
        dataStoreFile = dataAbsltLocation + "1_MCF.txt"
    elif formulation == "3":
        dataStoreFile = dataAbsltLocation + "1_STEINER.txt"
    elif formulation == "4":
        dataStoreFile = dataAbsltLocation + "1_NS.txt"

    if (UseLocalBranch == "0"):
        AveTime = 0.0
        f1 = open(str(dataStoreFile), "r")
        lines = f1.readlines()
        for i in range(1, int(graph_number) + 1):
            # begin collect data
            str_list = lines[-i].split()
            AveTime += (float(str_list[1]) / float(graph_number))
        f1.close()

        f1 = open(str(dataStoreFile), "a")
        _str = "Ave Time: " + str(round(AveTime, 3))
        print(_str)
        f1.write("-----------------HERE----------------------\n")
        f1.write(_str)
        f1.write("\n-----------------END----------------------\n\n")
        f1.close()

    else:
        Ave_TOT_TIME = 0.0
        Ave_LB_TIME = 0.0
        Ave_FINAL_SOLVE_TIME = 0.0
        Ave_LB_RUN_TIME = 0

        f1 = open(str(dataStoreFile), "r")
        lines = f1.readlines()
        for i in range(1, int(graph_number)+1):
            str_list = lines[-i].split()

            Ave_TOT_TIME += (float(str_list[1]) / float(graph_number))
            Ave_LB_TIME += (float(str_list[2]) / float(graph_number))
            Ave_FINAL_SOLVE_TIME += (float(str_list[3]) / float(graph_number))
            Ave_LB_RUN_TIME += (float(str_list[4]) / float(graph_number))

        f1.close()

        f1 = open(str(dataStoreFile), "a")
        _str1 = "Ave_TOT_TIME = " + str(round(Ave_TOT_TIME, 3))
        _str2 = "\nAve_LB_TIME = " + str(round(Ave_LB_TIME, 3))
        _str3 = "\nAve_FINAL_SOLVE_TIME = " + \
            str(round(Ave_FINAL_SOLVE_TIME, 3))
        _str4 = "\nAve_LB_RUN_TIME = " + str(round(Ave_LB_RUN_TIME))
        _str = _str1 + _str2 + _str3 + _str4

        print(_str)
        f1.write("-----------------HERE----------------------\n")
        f1.write(_str)
        f1.write("\n-----------------END----------------------\n\n")
        f1.close()
