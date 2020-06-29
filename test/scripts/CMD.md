`py smpTestNS.py  random_graph plan_random group_1 dataset1_1_1_1 1 4 3 0 1200 10`

`py smpTest.py  random_graph plan_random group_1 dataset1_1_3_2 33554431 4 3 0  0 3600 1200 0.2 1200 0.2`

`py smpParamAdjust.py  random_graph plan_random group_1 dataset1_1_3_2 1 4 3 0 1200 10`

`py smpTest_smallscale.py  random_graph plan_random group_1 tg 19 4 3 0  0 3600 1200 0.2 1200 0.2`

dense graph
py smpTest_smallscale.py  random_graph general_random group_1 n200_t30_p3_b50 11 4 3 0  0 3600 1200 0.2 1200 0.2

'py smpTest_smallcheck.py  random_graph plan_random group_1 tg 1_MCF 1_NS'

dense graph
'py smpTest_smallcheck.py  random_graph general_random group_1 n200_t30_p3_b50 1_MCF 1_NS'

- 25 1: 33554431
- 26th 1: 33554432