# %%
# Run this file at manytimes/ directory that contains manytimes/default/ directory.

import yaml
from cProfile import label
import os
import sys
from os import lseek, times
from tokenize import cookie_re
from xmlrpc.client import boolean
from numpy.core.shape_base import atleast_2d
from numpy.lib import type_check
from scipy import fftpack
import math
import numpy as np
import matplotlib.pyplot as plt
import csv
from scipy import signal
from python_modules import FFT_utility as fu

# %%
file_to_process = "beta_sum_corr.txt"
dir_index_start = int(1)
dir_index_end = int(50)
dirs_num = dir_index_end-dir_index_start+1
values_array = np.zeros(dirs_num)
print("Fetching values from "+str(dirs_num)+" directories...")

for dir_index in range(dir_index_start, dir_index_start+dirs_num):
    with open("./"+str(dir_index)+"/"+file_to_process) as raw_data:
        str_raw_data = raw_data.read()
    # correlate dir_index and values_array's index
    values_array_index = dir_index - dir_index_start
    values_array[values_array_index] = float(str_raw_data)

# %%
with open("./default/config.yaml", "r") as yml:
    config = yaml.safe_load(yml)
defects_rate = config["defects_rate"]
beta_mean = np.mean(values_array)
beta_err = np.std(values_array, ddof=1)/np.sqrt(dirs_num)
with open("beta_stat.csv", "w") as b:
    writer = csv.writer(b, delimiter=",")
    writer.writerow([defects_rate, beta_mean, beta_err])
print(defects_rate, beta_mean, beta_err)
# %%
