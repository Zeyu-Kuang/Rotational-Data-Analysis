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


## Motivation
The putative Lorentz violation would have several frequency components. 
The smallest one is the annual frequency of the earth.
This means that we need at least one year of data in order to decouple all the frequency components.
Since we are taking data at the frequency of 1 HZ on a rotating table, we would end up with billions of data.
We use linear fitting method to extract the frequency components. However, traditional fitting method will easily run out of the memory of any computer.
Motivated by this, we designed the demodulated least squre (DLS) fitting method to cope with big data.
DLS method firstly demodulates the large frequency components. Then step by step it goes to the slowest (such as annual frequency of the earth).
In this way, we can relieve the computational pressure without losing too much information (compared to the averaging).

## Summary
The program reduced the dataset using demodulated least square method. The result constraints the SME coefficient c_Q^T to 〖10〗^(-17).


## Explanation of the program

The 'Rotational Data Analysis' program is designed to analysis the measured data from the rotational experiment.

The 'Rotational Data Analysis' program is has two important features: calculate the SME coefficients in two methods: DLS and OLS, as well as the optimization of them. We explain the usage of important folders here.

\paragraph{./} We start by explaining the main folder\footnote{`./' means `Rotational Data Analysis/'}. It uses Data/ to store the raw data and processed data, as well as the result. It has Images/ to store generated image files. matlabProgram/ contains the same program, but written in MATLAB by Paul Louis Stanwix. tools/ is the local library used to store the shared functions of the program scripts. 

The rest of the files in \textbf{./} are .ipynb scripts. It means they are jupyter notebook scripts in python language. These scripts is for the main analysis. ReadingRotationData.ipynb is to read and visualize data in each runs. DLS-version2.ipynb uses DLS method to analyse data. OLS\textunderscore uses OLS method to slice and analyse the data. Both the optimization scripts are also in the same folder. These optimization scripts are used to optimize the parameters in DLS and OLS, such as subset rotational number, $N_r$, in DLS method. Plotting\textunderscore images.ipynb is used to generate quality figures for both presentation and report.
