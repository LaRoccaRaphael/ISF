#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import numpy as np
from scipy import interpolate
from scipy.signal import find_peaks
import os

def pick_picking(mzs_n,ints_n,distance,max_peaks):
    
    peaks, _ = find_peaks(ints_n,distance=distance)

    instensity_mean = ints_n[peaks]
    mz_peaks = mzs_n[peaks]

    peaks = peaks[np.argsort(ints_n[peaks])[::-1][:max_peaks]]
    selected_mzs = mzs_n[peaks]
    selected_intensity_arr = ints_n[peaks]

    new_mzs = selected_mzs[np.argsort(selected_mzs)]
    new_intensity_arr = selected_intensity_arr[np.argsort(selected_mzs)]
    
    return(new_mzs)


def list_files_with_suffix(directory, suffix):
    files = []
    for filename in os.listdir(directory):
        if filename.endswith(suffix):
            files.append(filename)
    return files



my_parser = argparse.ArgumentParser(allow_abbrev=False)
my_parser.add_argument('-i','--input', action='store', type=str, required=True,help='folder path containing mean spectrum files')
my_parser.add_argument('-t','--tolerance', action='store', type=float, required=True,help='Peak picking tolerance Da')
my_parser.add_argument('-o1','--output1', action='store', type=str, required=True,help='ouput file name')

args = my_parser.parse_args()
output_file = args.output1
tolerance = args.tolerance*10000


meanspectrum_f = list_files_with_suffix(args.input, 'meansp.npy')


mzs_n = np.arange(150, 1600, 0.0001)
all_mean = np.zeros(len(mzs_n))
    
for f in meanspectrum_f:
    ms = np.load(args.input+f)
    mzs = ms[0,]
    ints = ms[1,]
        
    ints = ints/np.max(ints)
    f = interpolate.interp1d(mzs , ints,fill_value="extrapolate")
    ints_n = f(mzs_n)
    ints_n = ints_n/np.sum(ints_n[ints_n>0])
    all_mean = all_mean + ints_n
        
ints_n = all_mean


new_mzs = pick_picking(mzs_n,ints_n,tolerance,10000)

np.savetxt(output_file,new_mzs)