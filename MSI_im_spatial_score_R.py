import numpy as np
import os
import argparse

my_parser = argparse.ArgumentParser(allow_abbrev=False)
my_parser.add_argument('-i','--input', action='store', type=str, required=True,help='input file msi')
my_parser.add_argument('-n','--noise', action='store', type=float, required=True,help='Noise intensity value')

args = my_parser.parse_args()

msi_file = args.input
noise = args.noise


d = np.load(msi_file+"_crd.npy" )
np.savetxt(msi_file+"_crd.csv",d,delimiter="," )
    
d = np.load(msi_file+"_msi.npy" )
d = d -noise 
d[d < 0] = 0

np.savetxt(msi_file+"_msi.csv",d,delimiter="," )

os.system("Rscript MSI_im_spatial_score.R " + msi_file)  


