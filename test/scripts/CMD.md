`py smpTestNS.py  random_graph plan_random group_1 dataset1_1_1_1 1 4 3 0 1200 10`

`py smpTest.py  random_graph plan_random group_1 dataset1_1_3_2 33554431 4 3 0  0 3600 1200 0.2 1200 0.2`

`py smpParamAdjust.py  random_graph plan_random group_1 dataset1_1_3_2 1 4 3 0 1200 10`

LB_MaxRestart, LB_MaxIter, Rmin, Rmax, BCSolNum, BCTime
`python smpTest_smallscale.py  random_graph plan_random group_1 tg 19 2 0 0 0 3 4 10 30 10 3 0 3600 1200 0.2 1200 0.2 1 1 1`

dense graph
py smpTest_smallscale.py  random_graph general_random group_1 tg 10 4 3 0 0 3 4 10 30 10 3 0 3600 1200 0.2 1200 0.2 1 1 1

'py smpTest_smallcheck.py  random_graph plan_random group_1 dataset1_1_t_g 0MCFres 0NSres'

dense graph
'py smpTest_smallcheck.py  random_graph general_random group_1 n200_t30_p3_b50 1_MCF 1_NS'

测试参数：
py testindexstrongcut.py random_graph plan_random group_1 tg 19 4 3 0 0 3 4 10 30 10 3 0 3600 1200 0.2 1200 0.2 0 0 0

py testindexusercut.py random_graph plan_random group_1 dataset1_1_3_2 50 4 3 0 0 3 4 10 30 10 3 0 3600 1200 0.2 1200 0.2 0 0 1

自动跑文件夹：
python autorun.py  random_graph general_random group_1 n1500_t30_p3_b0015_v08 20 4 1 0 1 3 4 10 30 10 3 0 3600 1200 0.2 1200 0.2 0 0 1

自动跑LB
python autorunLB.py  random_graph general_random group_1 n1500_t30_p3_b0015_v08 20 4 1 0 0 3 4 10 30 10 3 0 3600 1200 0.2 1200 0.2 1 0 1

跑Dimacs测试
run classic steiner
python runclassicsteiner.py  random_graph classic_steiner group_1 Copenhagen14 20 4 1 0 0 3 10 10 30 10 20 1 3600 1200 0.2 1200 0.2 1 0 1

跑Random测试
python runclassicsteiner.py  random_graph general_random group_1 Copenhagen14 20 4 1 0 0 3 10 10 30 10 20 1 3600 1200 0.2 1200 0.2 1 0 1

跑Plangraph测试
python runPlanGraph.py  random_graph plan_random group_1 Copenhagen14 20 2 0 0 0 3 10 10 30 10 20 1 3600 1200 0.2 1200 0.2 0 0 1

跑terminal
python runTerminal.py  random_graph general_random group_1 0-terminal 3 4 3 0 1 3 10 10 30 10 20 1 3600 1200 0.2 1200 0.2 0 0 1

跑Complex
python runComplex.py random_graph hardInstance group_1 auto 20 4 1 0 0 3 10 10 30 10 20 1 3600 1200 0.2 1200 0.2 0 0 1

跑cor
python runCoefficient.py  random_graph cor group_1 coefficients 20 4 1 0 0 3 10 10 30 10 20 1 3600 1200 0.0 1200 0.0 0 0 1