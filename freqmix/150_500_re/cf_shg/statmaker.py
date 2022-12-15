# %%
#coding: UTF-8
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

#%%
beta_mori=np.array([[0, 2.574,0],[0.01, 2.396,0]])
mori_0=2.574
mori_1=2.396
beta_shg=np.zeros((4,3))

file_to_process = "beta50_01percent.csv"
defects_rate=0.001
data=np.loadtxt(file_to_process,delimiter=',')
beta_mean = np.mean(data)
beta_err = np.std(data, ddof=1)/np.sqrt(len(data))
print(defects_rate, beta_mean, beta_err)
beta_shg[0]=defects_rate, beta_mean, beta_err
file_to_process = "beta50_02percent.csv"
defects_rate=0.002
data=np.loadtxt(file_to_process,delimiter=',')
beta_mean = np.mean(data)
beta_err = np.std(data, ddof=1)/np.sqrt(len(data))
print(defects_rate, beta_mean, beta_err)
beta_shg[1]=defects_rate, beta_mean, beta_err
file_to_process = "beta50_05percent.csv"
defects_rate=0.005
data=np.loadtxt(file_to_process,delimiter=',')
beta_mean = np.mean(data)
beta_err = np.std(data, ddof=1)/np.sqrt(len(data))
print(defects_rate, beta_mean, beta_err)
beta_shg[2]=defects_rate, beta_mean, beta_err
file_to_process = "beta50_07percent.csv"
defects_rate=0.007
data=np.loadtxt(file_to_process,delimiter=',')
beta_mean = np.mean(data)
beta_err = np.std(data, ddof=1)/np.sqrt(len(data))
print(defects_rate, beta_mean, beta_err)
beta_shg[3]=defects_rate, beta_mean, beta_err
beta_shg[3]=0, 2.56, 0


beta_mix=np.zeros((4,3))
file_to_process = "../0vac/manytimes/beta_stat.csv"
defects_rate=0.00
data=np.loadtxt(file_to_process,delimiter=',')
beta_mix[0]=data
file_to_process = "../03vac/manytimes/beta_stat.csv"
defects_rate=0.003
data=np.loadtxt(file_to_process,delimiter=',')
beta_mix[1]=data
file_to_process = "../05vac/manytimes/beta_stat.csv"
defects_rate=0.05
data=np.loadtxt(file_to_process,delimiter=',')
beta_mix[2]=data
file_to_process = "../07vac/manytimes/beta_stat.csv"
defects_rate=0.007
data=np.loadtxt(file_to_process,delimiter=',')
beta_mix[3]=data

fig,ax = plt.subplots()
ax.set_title("β vs monovacancy density")
ax.errorbar(beta_mix[:,0],beta_mix[:,1], yerr = beta_mix[:,2], capsize=5, fmt='o', markersize=12, label="mixing(This study)")
ax.errorbar(beta_shg[:,0],beta_shg[:,1], yerr = beta_shg[:,2], capsize=5, fmt='o', markersize=12, label="SHG(This study)")
ax.plot(beta_mori[:,0],beta_mori[:,1],marker="o", markersize=12, linestyle="none", label="SHG[5]")
ax.legend()
ax.set_xlabel("Vacancy density [-]")
ax.set_ylabel("β")
plt.savefig("vacrate_vs_beta_2ways.svg")
plt.show()


# %%
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
