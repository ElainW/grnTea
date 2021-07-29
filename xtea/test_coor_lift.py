#!/bin/python3
# 2021/07/27
##@@author: Yilan (Elain) Wang, Harvard University
##@@contact: yilanwang@g.harvard.edu
'''
1. testing coor_lift.py
'''

from coor_lift import *
sf_out = "/n/scratch3/users/y/yw222/xTEA_ML/test_ref/I0231_1_clip2/liftover_test.txt"
ref_bed = "/n/data1/bch/genetics/lee/elain/xTEA_benchmarking/alt_reference/gold_std_hg38_v4_int.temp"
coor_lift = Coor_Lift(sf_out, sf_out + ".lifted", ref_bed)
coor_lift.run_coor_lift("ctrl")