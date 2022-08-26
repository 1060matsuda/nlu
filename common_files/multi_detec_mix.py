import yaml
from cProfile import label
import os
import sys
from os import times
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

HPC_OR_LOCAL = "LOCAL"
imputFileName = "outp_1.csv"
data_raw = fu.loadCsvOutput(imputFileName)
freq_f1 = 100  # GHz
"""Bigger frequency [GHz]"""
freq_f2 = 30
"""Smaller frequency [GHz]"""
freq_fsum = freq_f1+freq_f2
"""Frequency sum [GHz]"""
freq_fdif = freq_f1 - freq_f2
"""Frequency difference [GHz]"""
freq_fgcd = np.gcd(freq_f1, freq_f2)
"""Greatest common divisor of frequencies f1 and f2 [GHz]"""
source_amp_f1 = 1000/freq_f1/20  # A
"""Amplitude of f1 wave [Å]"""
source_amp_f2 = 1000/freq_f2/20  # A
"""Amplitude of f2 wave [Å]"""
#freq_index = 0
non_detec_cols = int(3)
# Nc doesnt have to be integer
Nc_fgcd = int(1)
Nc_f1 = Nc_fgcd*(freq_f1/freq_fgcd)
Nc_f2 = Nc_fgcd*(freq_f2/freq_fgcd)  # How many cycles to window
v_expt = 5500  # m/s
"""Experimental value of the longitudinal wave velocity [m/s]"""

T_f1 = 1000/freq_f1
"""wave f1's wave cycle [ps]"""
T_f2 = 1000/freq_f2
"""wave f1's wave cycle [ps]"""
T_fgcd = 1000/freq_fgcd
"""Least common multiple of the wave cycles T1 and T2 [ps]. Represents the wave cycle of f_gcd."""

data_raw_cols = data_raw.shape[1]
detecs_num = data_raw_cols - non_detec_cols
# In most cases, data_raw[:,0] is time, data_raw[:,1] is source, data_raw[:,<last column>] is backend.

time = data_raw[:, 0]
timestep = time[1] - time[0]
if (int(time[-1]/timestep)!=len(time)+1):
    print("Data is possibly broken. Check if the data lacks information at some timestep.")
    sys.exit()
x_source = data_raw[:, 1]
x_detecs_array = data_raw[:, 2:2+detecs_num]
x_end = data_raw[:, -1]
# In most cases, data_raw[:,0] is time, data_raw[:,1] is source, data_raw[:,<last column>] is backend.
a_f1_at_detecs = np.zeros(detecs_num)
"""[m]"""
a_f2_at_detecs = np.zeros(detecs_num)
"""[m]"""
a_fsum_at_detecs = np.zeros(detecs_num)
"""[m]"""
a_fdif_at_detecs = np.zeros(detecs_num)
"""[m]"""

# not needed for freqmix
a_2f1_at_detecs = np.zeros(detecs_num)
"""[Å]"""
a_3f1_at_detecs = np.zeros(detecs_num)
"""[Å]"""
a_4f1_at_detecs = np.zeros(detecs_num)
"""[Å]"""

beta_SHGs_at_detecs = np.zeros(detecs_num)
beta_SHGs_at_detecs_corrected = np.zeros(detecs_num)
beta_sums_at_detecs = np.zeros(detecs_num)
beta_sums_at_detecs_corrected = np.zeros(detecs_num)
beta_difs_at_detecs = np.zeros(detecs_num)
beta_difs_at_detecs_corrected = np.zeros(detecs_num)
beta_aves_at_detecs= np.zeros(detecs_num)
#betaSHGArrayCorrected = np.zeros((len(freqs_f_1), len(detecs)))

#wavelength_f1 = np.zeros(len(freqs_f_1))
#waveVelocityArray = np.zeros((len(freqs_f_1), len(detecs)))

Ns_f1 = T_f1/timestep  # Ns: Number of data points in one cycle
# N: Total Number of data points in thw windowed region
N_f1 = int(Nc_f1*Ns_f1)

Ns_f2 = T_f2/timestep  # Ns: Number of data points in one cycle
# N: Total Number of data points in thw windowed region
N_f2 = int(Nc_f2*Ns_f2)

Ns_fgcd = T_fgcd/timestep
N_fgcd = int(Nc_fgcd*Ns_fgcd)

if(N_f1 != N_f2 or N_f2 != N_fgcd or N_fgcd != N_f1):
    print("REVIEW YOUR INPUT. Program has gone into forced shutdown.")
    sys.exit()

u_detecs_array = x_detecs_array - x_detecs_array[0]

delta_x_source_detec0 = x_detecs_array[0, 0] - x_source[0]
delta_x_detecs = x_detecs_array[0, 1] - x_detecs_array[0, 0]
delta_t_detecs = delta_x_detecs/10*(10**3)/5500  # ps
trimmed_waves = np.zeros((detecs_num, N_fgcd))
windowed_waves = np.zeros((detecs_num, N_fgcd))
wave_velocity_at_detecs_array = np.zeros(detecs_num)
zerocross_timestep_at_detecs_array = np.zeros(detecs_num)
arrival_timestep_at_detecs_array = np.zeros(detecs_num)

for i in range(detecs_num):
    x_detec = x_detecs_array[:, i]
    u_detec = u_detecs_array[:, i]
    delta_x_source_to_thisdetec = x_detec[0] - x_source[0]
    zerocross_timestep, arrival_timestep, wave_velocity = fu.wave_arrival_zerocross(
        u_detec, delta_x_source_to_thisdetec, v_expt, timestep, T_fgcd)

    wave_velocity_at_detecs_array[i] = wave_velocity
    zerocross_timestep_at_detecs_array[i] = zerocross_timestep
    arrival_timestep_at_detecs_array[i] = arrival_timestep

wave_velocity = wave_velocity_at_detecs_array[-1]
del zerocross_timestep
del arrival_timestep

wavelength_f1 = wave_velocity * T_f1 * (10**(-12))
"""[m]"""
wavelength_f2 = wave_velocity * T_f2 * (10**(-12))
"""[m]"""
k_f1 = 2*np.pi/wavelength_f1
"""[1/m]"""
k_f1 = 2*np.pi/wavelength_f1
"""[1/m]"""

# フーリエ変換
for i in range(detecs_num):
    # ゼロクロス法の音速測定
    # Δx_source_and_detec[A]*100[pm/A] / v[pm/ps] = Δt_source_and_detec [ps]
    # この時刻から半周期たった時（first positive peak）を基準に，そこからゼロクロス法を行う
    x_detec = x_detecs_array[:, i]
    u_detec = u_detecs_array[:, i]
    delta_x_source_to_thisdetec = x_detec[0] - x_source[0]

    window_start_timestep = int(
        zerocross_timestep_at_detecs_array[i] - 2*T_f2/timestep)

    # ここまで音速測定

    trimmedWave = u_detec[window_start_timestep:window_start_timestep+N_fgcd]
    trimmedTime = time[window_start_timestep:window_start_timestep+N_fgcd]

    windowedWave = fu.window(trimmedWave, "hann")
    #paddedWave = fu.zeroPadding(windowedWave)
    #paddedWaveWithoutWindow = fu.zeroPadding(trimmedWave)

    ##SUPER IMPORTANT##
    waveToTransform = windowedWave
    ##CHOU DAIJI##

    """if (i==detecs_num-2):
        plt.plot(time, u_detec)
        plt.axvline(x=window_start_timestep*timestep)
        plt.axvline(x=(window_start_timestep+N_fgcd)*timestep)
        plt.show()"""

    # FFT. transformedArray: [0]=power, [1]=freq
    # FFTedData = fftWithWindow(trimmedWave, "hann") #window = "hann" or "hamming"
    #FFTedData = fftWithWindow(zeroPadding(trimmedWave), "hann")
    ffted_data = fu.FFTonly(waveToTransform, timestep)
    #FFTedData=FFTonly(window(trimmedWave, "hann"))
    abs_ffted_data = np.abs(ffted_data)

    index_fsum = fu.getIndexOfNearestValue(abs_ffted_data[1], freq_fsum*10**9)
    A_fsum = abs_ffted_data[0][index_fsum]
    a_fsum = A_fsum*2/int(len(waveToTransform))*10**-10
    a_fsum_at_detecs[i] = a_fsum

    index_fdif = fu.getIndexOfNearestValue(abs_ffted_data[1], freq_fdif*10**9)
    A_fdif = abs_ffted_data[0][index_fdif]
    a_fdif = A_fdif*2/int(len(waveToTransform))*10**-10
    a_fdif_at_detecs[i] = a_fdif

    index_f1 = fu.getIndexOfNearestValue(abs_ffted_data[1], freq_f1*10**9)
    A_f1 = abs_ffted_data[0][index_f1]
    a_f1 = A_f1*2/int(len(waveToTransform))*10**-10
    a_f1_at_detecs[i] = a_f1

    index_f2 = fu.getIndexOfNearestValue(abs_ffted_data[1], freq_f2*10**9)
    A_f2 = abs_ffted_data[0][index_f2]
    a_f2 = A_f2 * 2/int(len(waveToTransform))*10**-10
    a_f2_at_detecs[i] = a_f2

    beta_mix = fu.getBetaFreqMix(a_fsum, a_fdif, source_amp_f1*10**-10, source_amp_f2*10**-10,
                                 freq_f1*10**9, freq_f2*10**9, delta_x_source_to_thisdetec*10**-10, wave_velocity)
    beta_mix_corrected = fu.getBetaFreqMix(
        a_fsum, a_fdif, a_f1, a_f2, freq_f1*10**9, freq_f2*10**9, delta_x_source_to_thisdetec*10**-10, wave_velocity_at_detecs_array[i])

    beta_sum = beta_mix[0]
    beta_sum_corrected = beta_mix_corrected[0]
    beta_sums_at_detecs[i] = beta_sum
    beta_sums_at_detecs_corrected[i] = beta_sum_corrected

    beta_dif = beta_mix[1]
    beta_dif_corrected = beta_mix_corrected[1]
    beta_difs_at_detecs[i] = beta_dif
    beta_difs_at_detecs_corrected[i] = beta_dif_corrected

    # higher harmonics amplitude[arb]
    harmonicsIndex = fu.getIndexUpToSixthHarmonic(
        abs_ffted_data[1], freq_f1*10**9)
    A_f1 = abs_ffted_data[0][harmonicsIndex[0]]
    A_2f1 = abs_ffted_data[0][harmonicsIndex[1]]
    A_3f1 = abs_ffted_data[0][harmonicsIndex[2]]
    A_4f1 = abs_ffted_data[0][harmonicsIndex[3]]
    A_5f1 = abs_ffted_data[0][harmonicsIndex[4]]
    A_6f1 = abs_ffted_data[0][harmonicsIndex[5]]

    a_f1 = A_f1*2/int(len(waveToTransform))*10**-10
    a_f1_at_detecs[i] = a_f1

    a_2f1 = A_2f1*2/int(len(waveToTransform))*10**-10
    a_2f1_at_detecs[i] = a_2f1

    a_3f1 = A_3f1*2/int(len(waveToTransform))*10**-10
    a_3f1_at_detecs[i] = a_3f1

    a_4f1 = A_4f1*2/int(len(waveToTransform))*10**-10
    a_4f1_at_detecs[i] = a_4f1

    a_f1_source = source_amp_f1*10**-10
    a_f2_source = source_amp_f2*10**-10

    beta_shg = fu.getBetaSHG(a_f1_source, a_2f1, wavelength_f1,
                             delta_x_source_to_thisdetec*(10**-10))
    beta_shg_corrected = fu.getBetaSHG(
        a_f1, a_2f1, wavelength_f1, delta_x_source_to_thisdetec*(10**-10))
    """with open("betaSHG.txt", "w") as f:
        f.write(str(betaSHG))
    """
    beta_SHGs_at_detecs[i] = beta_shg
    beta_SHGs_at_detecs_corrected[i] = beta_shg_corrected

    if i == detecs_num-3:
        fig, ax = plt.subplots()
        ax.plot(trimmedTime, waveToTransform)
        plt.savefig("wave_of_detector #"+str(i)+".png", format="png")
        

beta_aves_at_detecs=(beta_sums_at_detecs+ beta_difs_at_detecs)/2

fig, ax = plt.subplots()
ax.set_xlabel(r"propagation distance $x (\mathrm{\AA)}$ ")
ax.set_ylabel(r"wave velocity $v (m/s)$")
ax.plot(x_detecs_array[0], wave_velocity_at_detecs_array,
        marker="o", linestyle="none")
plt.savefig("x_vs_velocity.png")

fig, ax = plt.subplots()
title = "x vs f1 amplitude"
plt.plot(x_detecs_array[0, :], a_f1_at_detecs, marker="o", linestyle="none")
plt.title(title)
plt.savefig(title+".png")

fig, ax = plt.subplots()
title="x vs f2 amplitude"
plt.plot(x_detecs_array[0, :], a_f2_at_detecs, marker="o", linestyle="none")
plt.title(title)
plt.savefig(title+".png")

fig, ax = plt.subplots()
title = "x vs fsum amplitude"
plt.plot(x_detecs_array[0, :], a_fsum_at_detecs, marker="o", linestyle="none")
plt.title("x vs fsum amplitude")
plt.savefig(title+".png")

fig, ax = plt.subplots()
title="x vs fdif amplitude"
plt.plot(x_detecs_array[0, :], a_fdif_at_detecs, marker="o", linestyle="none")
plt.title(title)
plt.savefig(title+".png")


fig, ax = plt.subplots()
ax.set_title("$a/x$")
ax.plot(x_detecs_array[0, :], a_fsum_at_detecs /
        x_detecs_array[0, :], marker="o", linestyle="none", label=r"$a_{sum}/x$")
ax.plot(x_detecs_array[0, :], a_fdif_at_detecs /
        x_detecs_array[0, :], marker="o", linestyle="none", label=r"$a_{dif}/x$")
ax.legend()
plt.savefig("x_vs_a.png")


# slope method
y_sum_source=  a_fsum_at_detecs/a_f1_source/a_f2_source
fig,ax=plt.subplots()
ax.plot(x_detecs_array[0, :], y_sum_source[:], marker="o", linestyle="none", label="source")

y_sum_corrected = a_fsum_at_detecs/a_f1_at_detecs/a_f2_at_detecs
ax.plot(x_detecs_array[0, :], y_sum_corrected[:], marker="o", linestyle="none", label="corrected")
title="y_sum"
ax.set_title(title)
ax.legend()
plt.savefig(title+".png")

fig,ax=plt.subplots()
y_dif_source=  a_fdif_at_detecs/a_f1_source/a_f2_source
ax.plot(x_detecs_array[0, :], y_dif_source[:], marker="o", linestyle="none", label="source")

y_dif_corrected = a_fdif_at_detecs/a_f1_at_detecs/a_f2_at_detecs
ax.plot(x_detecs_array[0, :], y_dif_corrected[:], marker="o", linestyle="none", label="corrected")
title="y_dif"
ax.set_title(title)
ax.legend()
plt.savefig(title+".png")
# detector 0からM-1 が外れ値をとっている場合がある．その場合N以降でslopeを取ったほうが良い．
# 外れていない場合は M=0
