# %%
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
with open("target_folders_plotter.csv", "r") as f:
    reader = csv.reader(f)
    dirs = [row for row in reader]
dirs = dirs[0]
print("Directory search completed")
print(dirs)

# %%
dirs_num = len(dirs)
rates_array = np.zeros(dirs_num)
betas_array = np.zeros(dirs_num)
errs_array = np.zeros(dirs_num)
print("Data input...")
for dir_name in dirs:
    print("from"+str(dir_name))
    index = dirs.index(dir_name)
    target_file = dir_name+"manytimes/beta_stat.csv"
    with open(target_file, "r") as f:
        reader = csv.reader(f)
        values = [row for row in reader]
    values = values[0]
    rates_array[index] = float(values[0])*100
    betas_array[index] = float(values[1])
    errs_array[index] = float(values[2])

print("plotting...")
fig, ax = plt.subplots()
ax.errorbar(rates_array, betas_array, yerr=errs_array,
            linestyle="none", marker="o", capsize=3)
ax.set_xlabel(r"Density [%]")
ax.set_ylabel(r"Nonlinearity parameter $\beta$")
ax.set_title(r"Vacancy density VS $\beta$")
plt.savefig("vac_rate_vs_beta.svg")
