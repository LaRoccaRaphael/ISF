#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import numpy as np
from numpy import genfromtxt
from pyimzml.ImzMLParser import ImzMLParser

def create_mean_spectrum(msi_file,min_mz,max_mz):
    p = ImzMLParser(msi_file)
    mean_spectrum_mz = np.arange(0, max_mz, 0.0001)
    mean_spectrum_int = np.zeros(np.size(mean_spectrum_mz,0))
    i = 0 
    for idx, (x,y,z) in enumerate(p.coordinates):
        mzs, intensities = p.getspectrum(idx)
        intensities_norm = intensities/np.sum(intensities)
        mzs_ind = np.rint(mzs*10000).astype('int')
        mean_spectrum_int[mzs_ind] += intensities_norm
        
    mzs = np.where(mean_spectrum_int >0)[0]/10000
    intensities = mean_spectrum_int[np.where(mean_spectrum_int >0)[0]]
    return(mzs,intensities)


my_parser = argparse.ArgumentParser(allow_abbrev=False)
my_parser.add_argument('-i','--input', action='store', type=str, required=True,help='imzML file path')
my_parser.add_argument('-o1','--output1', action='store', type=str, required=True,help='output file path')


args = my_parser.parse_args()

msi_file = args.input
output_file = args.output1

min_mz = 0
max_mz = 5000

mzs, ints = create_mean_spectrum(msi_file ,min_mz,max_mz)
t =  np.vstack((mzs, ints))
np.save(output_file+"_meansp", t)