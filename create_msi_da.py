#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import numpy as np
from numpy import genfromtxt
from pyimzml.ImzMLParser import ImzMLParser
from scipy.signal import find_peaks
import scipy.signal as signal

def find_msi_size(p):
    i = 0 
    for idx, (x,y,z) in enumerate(p.coordinates):
        i += 1 
    return(i)

def check_mass_database(experimental_mass, database_mass, tolerance):
    if database_mass != 0:
        return abs((experimental_mass - database_mass)) <= tolerance
    
    
def binarySearch_tol(arr, l, r, x, tolerance): 
  
    while l <= r: 
  
        mid = l + (r - l)//2; 
          
        if check_mass_database(x,arr[mid],tolerance): 
            itpos = mid +1
            itneg = mid -1
            index = []
            index.append(mid)
            if( itpos < len(arr)):
                while check_mass_database(x,arr[itpos],tolerance) and itpos < len(arr):
                    index.append(itpos)
                    itpos += 1 
            if( itneg > 0): 
                while check_mass_database(x,arr[itneg],tolerance) and itneg > 0:
                    index.append(itneg)
                    itneg -= 1     
            return index 
  
        elif arr[mid] < x: 
            l = mid + 1
 
        else: 
            r = mid - 1
      
    return -1

def align_peaks(peaks_mz,intensities,database_exactmass, tolerance): 
    pixel_new_vec = np.zeros(np.size(database_exactmass,0))
    for i in range(0,np.size(peaks_mz,0)):
        exp_peak = peaks_mz[i]
        
        vec_ind = binarySearch_tol(np.append(database_exactmass,np.max(database_exactmass)+1), 0, len(database_exactmass)-1, exp_peak,tolerance)
        if vec_ind != -1:
            if len(vec_ind) > 1:
                for j in range(0,np.size(vec_ind,0)):
                    if intensities[i] > pixel_new_vec[vec_ind[j]]:
                        pixel_new_vec[vec_ind[j]] = intensities[i]
                        
            if len(vec_ind) == 1:
                if intensities[i] > pixel_new_vec[vec_ind]:
                    pixel_new_vec[vec_ind] = intensities[i]
                
    return(pixel_new_vec)


my_parser = argparse.ArgumentParser(allow_abbrev=False)
my_parser.add_argument('-i','--input', action='store', type=str, required=True,help='path to msi')
my_parser.add_argument('-i2','--input2', action='store', type=str, required=True,help='path to file containing m/zs')
my_parser.add_argument('-o1','--output1', action='store', type=str, required=True,help='output path')

my_parser.add_argument('-tol','--tolerance', action='store', type=float, required=True,help='m/z tolerance in da for the identification')

args = my_parser.parse_args()

msi = args.input
peaks_file = args.input2
output_file = args.output1

tolerance = args.tolerance



msi_mz = genfromtxt(peaks_file, delimiter=' ')
database_exactmass = np.sort(msi_mz)

p = ImzMLParser(msi)

pixel_number = find_msi_size(p)
msi = np.zeros((np.size(database_exactmass),pixel_number))
crd = np.zeros((pixel_number,2),dtype=int)
tic = np.zeros((pixel_number,2),dtype=int)

i = 0 
for idx, (x,y,z) in enumerate(p.coordinates):
    mzs, intensities = p.getspectrum(idx)
    intensities = np.array(intensities)
    peaks_mz = mzs
    vec_int = intensities
    tic[i,0] = np.sum(intensities)
    tic[i,1] = np.sum(vec_int)
    pixel_new_vec = align_peaks(peaks_mz,vec_int,database_exactmass, tolerance)
    msi[:,i] = pixel_new_vec
    crd[i,0] = x
    crd[i,1] = y
    i += 1


np.save(output_file +'_msi', msi)
np.save(output_file +'_crd', crd)
np.save(output_file +'_tic', tic) 