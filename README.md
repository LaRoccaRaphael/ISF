# Image Structure Filtering (ISF) method

The method automatically compares the ion images from MSI belonging to two conditions (control and interaction), and identifies the specific ions for each condition. One unique aspect of the ISF is that it considers the difference in the spatial distribution of the same ion in different MSI conditions by assessing the topological coherence of each ion image.

The ISF method involves 5 steps:

1. Generate the average spectrum of all MSI or a subset of MSI.
2. Combine all MSI average spectra, detect the peaks, and find all m/z local maxima.
3. Generate MSI datacube by creating an ion image for each detected peak in the average spectrum.
4. Assign a score based on the topological coherence of each ion image in each MSI.
5. Compare the scores of the ion images of the same ion between different conditions.


This repository contains the scripts for all the steps of the method developed in the paper. The MSI is recalibrated independently by (https://github.com/LaRoccaRaphael/MSI_recalibration.git)



## ISF: how to run the code? 

### Prerequisites 
The code is written in <code>python version 3.6.10</code> and runs with the following libraries: 

```
numpy==1.18.1 
argparse==1.1
pyimzml==1.5.1
scipy== 1.5.4
pandas==1.1.5
```

One script is written in <code>R version 4.2.1</code> with the following library: 
    
```
mmand== 1.6.3
```

### Command lines

Here we explain how to use the scripts to execute each of the first 4 steps of the ISF method.

(1)
```
python ./mean_spectrum_msi.py -i /path/to/the/imzML/file.imzML -o1 /path/to/the/output/file/MSI_name
```

(2)
```
python ./peak_picking.py -i /path/to/dir/of/mean_spectra/ -t 0.01 -o1 /path/to/the/list/of/masses/file.txt
```
The value of the following argument <code>-t</code> corresponds to the window size used during peak picking.

(3)
```
python ./create_msi_da.py -i /path/to/the/imzML/file.imzML -i2 /path/to/the/list/of/masses/file.txt  -tol 0.01 -o1 /path/with/new/MSI/name
```
The value of the following argument <code>-tol</code> corresponds to the m/z window around the m/z from the peak list from which we report the maximum intensity value detected in each pixel of an MSI. 

(4)
```
python ./MSI_im_spatial_score_R.py -i /path/with/new/MSI/name -n 10000
```

The value of the following argument <code>-n</code> corresponds to an absolute intensity cut-off below which we consider every value as 0. 

Note that you will need to replace */path/to/the/imzML/file.imzML*  by the path of the corresponding imzML file subjects to analysis; you will need to replace  */path/to/the/output/file/MSI_name* by the path to the output mean spectrum, *_mean_sp.npy* will be added automatically; you will need to replace */path/to/the/list/of/masses/file.txt* by the path of the file containing the values of the mz of the detected peaks; you will need to replace  */path/with/new/MSI/name* by the path to new MSI, the */name* should be without extension since they will be added automatically. The .imzML file and its corresponding .ibd file should be located in the same folder.
    
### Jupyter notebook Data analysis (step 5)
    
The notebook named <code>Data_analysis.ipynb</code> contains an example of how to load the processed MSI files and to compute the structure fold change as described in the paper.




## How to cite us?
Please cite our publication by using the following bibtex entry:



## License

Apache v2.0
See the [LICENSE](LICENSE) file for details.
