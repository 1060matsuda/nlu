from cProfile import label
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


def getDownwardZeroCrossIndex(vector1d):
    """searches for the first zerocross point (+ -> -) 
    from the leftend of the waveform."""
    downCount = 0
    searchIndex = 1
    while True:
        searchIndex = searchIndex + 1
        downCount = 0
        if vector1d[searchIndex] < 0 and vector1d[searchIndex-1] > 0:
            for i in range(10):
                if vector1d[searchIndex-5+i] < vector1d[searchIndex-5+i+1]:
                    downCount = 0
                else:
                    downCount = downCount+1
            if downCount == 10:
                # print(searchIndex)
                return searchIndex


def getUpwardZeroCrossIndex(vector1d):
    """searches for the first zerocross
    point (+ -> -) from the leftend of the waveform."""
    upCount = 0
    searchIndex = 1
    while True:
        upCount = 0
        searchIndex = searchIndex + 1
        if vector1d[searchIndex] > 0 and vector1d[searchIndex-1] < 0:
            for i in range(10):
                if vector1d[searchIndex-5+i] > vector1d[searchIndex-5+i+1]:
                    upCount = 0
                else:
                    upCount = upCount+1
            if upCount == 10:
                # print(searchIndex)
                return searchIndex


def getDownwardZeroCrossIndexFromArbitraryPoint(vector1d, startIndex):
    """Searches for an downward zerocross point 
    from the specified index of the vector."""
    downCount = 0
    searchIndex = startIndex
    while True:
        searchIndex = searchIndex + 1
        if vector1d[searchIndex] < 0 and vector1d[searchIndex-1] > 0:
            downCount = 0
            for i in range(10):
                if vector1d[searchIndex-5+i] < vector1d[searchIndex-5+i+1]:
                    downCount = 0
                else:
                    downCount = downCount+1
            if downCount == 10:
                # print(searchIndex)
                return searchIndex


def getUpwardZeroCrossIndexFromArbitraryPoint(vector1d, startIndex):
    """Searches for an upward zerocross point 
    from the specified index of the vector."""
    upCount = 0
    searchIndex = startIndex

    while True:
        upCount = 0
        searchIndex = searchIndex + 1
        if vector1d[searchIndex] > 0 and vector1d[searchIndex-1] < 0:
            for i in range(10):
                if vector1d[searchIndex-5+i] > vector1d[searchIndex-5+i+1]:
                    upCount = 0
                else:
                    upCount = upCount+1
                if upCount == 10:
                    # print(searchIndex)
                    return searchIndex


# Reads the csv file where time series data of detector coorinate is saved.
def loadCsvOutput(csvData):
    """Reads the output csv file and covert it into ndarray."""
    ndarrayData = np.loadtxt(csvData, delimiter=",")

    # Several timesteps are recorded twice in the raw csv file.
    # Below is the script which deletes the double-recorded timesteps.
    timeAtTheRowAbove = -100
    delRows = []
    for i in range(len(ndarrayData)):
        if timeAtTheRowAbove == ndarrayData[i, 0]:
            delRows.append(i)
        timeAtTheRowAbove = ndarrayData[i, 0]
    ndarrayData = np.delete(ndarrayData, delRows, 0)

    return ndarrayData

# Withdraws the time series of the detector's x coorinate


def getDetector(array):
    return array[:, 2]

# Withdraws the time series of the source's x coorinate


def getSource(array):
    return array[:, 1]

# trims the wave within the specified range


def trim(wave, trimStartIndex, trimRange):
    """trims the wave from the specified index, for the specified range"""
    trimmedArray = wave[trimStartIndex:trimStartIndex+trimRange]
    return trimmedArray

# trims the wave within the specified range and then offsets the wave so that the normal position can be x=0.


def trimAndOffset(wave, trimStartIndex, trimRange):
    """trims the wave within the specified range and then 
    offsets the wave so that the normal position can be x=0."""
    trimmedArray = wave[trimStartIndex:trimStartIndex+trimRange]
    offsettedArray = trimmedArray - wave[0]
    return offsettedArray

# performs fft. output[0] = power, output[1] = freq. both output are recognized as complex.


def fftWithWindow(FFTData, dtps, hannORhamming: str = "hann"):
    """Performs fft. Output[0] = power, output[1] = freq. 
    Both outputs are recognized as complex number.
    Timestep dt should be [ps] unit."""
    dataPoints = len(FFTData)
    if hannORhamming == "hann":
        windowFunction = signal.hann(dataPoints)
    elif hannORhamming == "hamming":
        windowFunction = signal.hamming(dataPoints)
    else:
        raise Exception(
            "Window is not/wrongly specified. Either hann / hamming is right.")
    acf = 1/(sum(windowFunction)/dataPoints)
    # print("acf")
    # print(acf)
    waveToTransform = acf*windowFunction*FFTData
    FFT_power = np.fft.fft(waveToTransform, n=None, norm=None)
    FFT_freq = np.fft.fftfreq(dataPoints, d=dtps*(10**-12))
    return np.stack([FFT_power, FFT_freq])


def window(data, hannORhamming: str = "hann"):
    dataPoints = len(data)
    if hannORhamming == "hann":
        windowFunction = signal.hann(dataPoints)
    elif hannORhamming == "hamming":
        windowFunction = signal.hamming(dataPoints)
    else:
        raise Exception(
            "Window is not/wrongly specified. Either hann / hamming is right.")
    acf = 1/(sum(windowFunction)/dataPoints)
    waveToTransform = acf*windowFunction*data
    return waveToTransform


def FFTonly(data, dtps):
    """Performs FFT without hann/hamming window."""
    FFT_power = np.fft.fft(data, n=None, norm=None)
    FFT_freq = np.fft.fftfreq(len(data), d=dtps*(10**-12))
    return np.stack([FFT_power, FFT_freq])

# calculates beta.


def getBetaSHG(a_f, a_2f, lambda_f, x_detec):
    """Calculates beta with SHG method. 
    Requires a_f, a_2f, lambda_f, x_D (distance between the source and the detector)."""
    beta = 8*a_2f*lambda_f*lambda_f/x_detec/a_f/a_f/np.pi/np.pi/2/2
    return beta

# attaches zeros to the raw wave data.


def zeroPadding(data):
    """Attaches zeros to the raw wave data."""
    zeros = np.zeros(len(data))
    buffer = 0
    for i in range(50):
        if len(data) < 2**i:
            buffer = 2**(i+4) - len(data)
            break
    zeros = np.zeros(int(buffer/2))
    paddedData = np.hstack((zeros, data, zeros))
    acf = (sum(np.abs(data)) / len(data)) / \
        (sum(np.abs(paddedData)) / len(paddedData))
    return acf*paddedData

# searches the 1d array data for the input value,
# and returns the index of the array where the nearest value of the input value is contained.


def getIndexOfNearestValue(data, value):
    """Returns the index i such that data[i] is the nearest to the input value.
    The input data should be one-dimensional."""
    index = np.argmin(np.abs(np.array(data) - value))
    return index

# searches the FFT-performed 1d array of frequencies, and returns the index of 1st - 6th harmonics.


def getIndexUpToSixthHarmonic(data, fundamental_frequency):
    """Returns the indices i1, i2, ... , i6 such that data[i1] is the nearest to fundamental_frequency, 
    data[i2] is the nearest to fundamental_frequency*2, ... , 
    data[i6] is the nearest to fundamental_frequency*6"""
    index1 = getIndexOfNearestValue(data, fundamental_frequency)
    index2 = getIndexOfNearestValue(data, fundamental_frequency*2)
    index3 = getIndexOfNearestValue(data, fundamental_frequency*3)
    index4 = getIndexOfNearestValue(data, fundamental_frequency*4)
    index5 = getIndexOfNearestValue(data, fundamental_frequency*5)
    index6 = getIndexOfNearestValue(data, fundamental_frequency*6)
    return np.array([index1, index2, index3, index4, index5, index6], dtype=np.int64)


def getBetaFreqMix(aSum, aDif, aF1, aF2, freq1, freq2, DeltaX, vel):
    """Calculates beta with freqmix method."""
    aMix = (aSum+aDif)/2
    Lambda1 = vel/freq1
    Lambda2 = vel/freq2
    K1 = 2*math.pi/Lambda1
    K2 = 2*math.pi/Lambda2
    return np.array([4*aDif/DeltaX/aF1/aF2/K1/K2, 4*aSum/DeltaX/aF1/aF2/K1/K2, 4*aMix/DeltaX/aF1/aF2/K1/K2])


def wave_arrival_zerocross(u_at_x, dx_wavesource_and_x, v_temp, timestep, T):
    """Returns the velocity calculated with zerocross method at x=x.
    You need to input displacement timeseries at x[A], 
    distance between the wavesource and x[A],
    temporary velocity (e.g. expt value)[m/s],
    timestep [ps/step], and the wave cycle T [ps].
    """
    dt_wavesource_and_x = dx_wavesource_and_x*100/v_temp  # ps
    zerocross_start_timestep = int(dt_wavesource_and_x/timestep)
    zerocross_timestep = getDownwardZeroCrossIndexFromArbitraryPoint(
        u_at_x, zerocross_start_timestep)
    arrival_timestep = int(zerocross_timestep - T/timestep/2)
    # v = dx [A] / (arrival_timestep[step]*timestep[ps/step])
    #   = dx/(arrival_timestep*timestep) [A/ps= 10^-10m/10^-12s = 100m/s]
    #   = dx/(arrival_timestep*timestep)*100 [m/s]
    velocity = dx_wavesource_and_x/(arrival_timestep*timestep)*100
    return zerocross_timestep, arrival_timestep, velocity