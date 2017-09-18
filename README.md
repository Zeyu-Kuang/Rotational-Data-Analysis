# Rotational Data Analysis
We are trying to find putative violation to Lorentz invariance in phonon section.
We build a rotating table with two phonon oscillators on top.
This program read the frequency difference between two phonon oscillators and determine the possible violation coefficients using two fitting methods.

## Background
Standard Model Extension (SME) was created to solve the incompatibility between our understanding of particle physics and gravitational physics. 
It allows potential violation of Lorentz invariance. 
We tested the Lorentz invariance of phonon using quartz oscillators. 
We generate precise frequencies from two quartz oscillators on a rotational table. 
The signal we measured is the frequency difference between the oscillator. 
Significant difference would result in a possible violation to isotropy of speed of the phonons.

Rotational Table: 

<img src="https://github.com/kzy2142/Rotational-Data-Analysis/blob/master/Images/rotation_experiment.png" width="200" height="150">


More information can be found in my summary: https://drive.google.com/open?id=0B0WjPWir-kYha2x6T3lpMGp5eFE


## Motivation
The putative Lorentz violation would have several frequency components. 
The smallest one is the annual frequency of the earth.
This means that we need at least one year of data in order to decouple all the frequency components.
Since we are taking data at the frequency of 1 HZ on a rotating table, we would end up with **billions of data**.
We use linear fitting method to extract the frequency components. However, traditional fitting method will easily run out of the memory of any computer.
Motivated by this, we designed the demodulated least squre (DLS) fitting method to cope with big data.
DLS method firstly demodulates the large frequency components. Then step by step it goes to the slowest (such as annual frequency of the earth).
In this way, we can relieve the computational pressure without losing too much information (compared to the averaging).

## Results
In this work, we analyzed the recorded data of 150 days using DLS method and restrict the some SME coefficients to the limit of 10^-17.
Further analysis of  data of a year length will generate more accurate result 

Optimization of DLS in run 8&9 and 11:

<img src="https://github.com/kzy2142/Rotational-Data-Analysis/blob/master/Images/Optimization.png" width="450" height="225">

We also studied the noise level in the data at different frequencies:

<img src="https://github.com/kzy2142/Rotational-Data-Analysis/blob/master/Images/compareRunsSubplots.png" width="600" height="300">

A potential source of fluctuation of the current experiment is the magnetic field in the laboratory. Further effort has been put on this area to monitor the magnetic field in the laboratory.

## Explanation of the program

The 'Rotational Data Analysis' program is designed to analysis the measured data from the rotational experiment.

The 'Rotational Data Analysis' program is has two important features: calculate the SME coefficients in two methods: DLS and OLS, as well as the optimization of them. We explain the usage of important folders here.

**./** We start by explaining the main folder('./' means 'Rotational Data Analysis/'). It uses Data/ to store the raw data and processed data, as well as the result. It has Images/ to store generated image files. matlabProgram/ contains the same program, but written in MATLAB by Paul Louis Stanwix. tools/ is the local library used to store the shared functions of the program scripts. 

The rest of the files in **./** are .ipynb scripts. It means they are jupyter notebook scripts in python language. These scripts is for the main analysis. ReadingRotationData.ipynb is to read and visualize data in each runs. DLS-version2.ipynb uses DLS method to analyse data. OLS_slicing uses OLS method to slice and analyse the data. Both the optimization scripts are also in the same folder. These optimization scripts are used to optimize the parameters in DLS and OLS, such as subset rotational number in DLS method. Plotting\textunderscore images.ipynb is used to generate quality figures for both presentation and report.

**./Data/** This folder contains raw data and processed data, as well as the result. The raw data is stored in the subfolder: original/. We record the pathname of the raw data of each run in a .txt file so that the program can read the raw data flexibly. For example, we record the pathname of run 8 in File_location_run_8.txt and use the program to read the file location of run 8, then read the raw data. Besides the raw data, **./Data/** folder also contains the processed data and the result. Since these data are produced by different .ipynb script, they are saved in different subfolders. Most subfolders have the same name as the script that generates this results. 

**./matlabProgram/** This folder contains the MATLAB program DLS1.m, used to calculate parameters of the first stage of DLS. Its result is stored in the subfolder Result/. We compare the result generated from MATLAB program and python program and latter has less numerical outliers.

**./tools/** This is the most important folder of this program because it contains all the necessary functions to perform the DLS and OLS methods. Inside this folder, dls.py and ols.py are two .py files. Since they are .py files, we can import them as library in python language. Apart from .py files, we also have .ipynb files in the same folder to explain the use of each function in the library. 
