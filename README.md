# Minimal Interspecies Interaction Adjustment (MIIA)
MIIA is a tool for predicting member-dependent interactions in microbial communities, which can help our deeper understanding of the organization principles in ecological communities.

The manuscript of MIIA has been submitted for publication and This repository contains the supplementary materials used to reproduce the simulation results submitted to Bioinformatics for review. 

## System requirements
Matlab R2015b+ in Windows, Mac OS, Linux

## Installation
1. Download
2. Unzip ./data/glvData.zip under ./data if you want to test ``glv`` datasets. 
3. Add path as follows:
```matlab
addpath(/path/to/miia/runMiia.m);
```

## Tutorial
Please describe the tutorial and how to get started. 

```matlab
run runMiia
```
You can test with a different data file by comment/uncomment the line 7 to 9.
```matlab
% dataSource = 'tutorial';
dataSource = 'friedman';  
% dataSource = 'glv';
```
You can reproduce the same results with the manuscript. Please refer to the below for each dataset.

## Datasets
1. Tutorial data
```matlab
dataSource = 'tutorial'
```
Figure 1.

2. Friedman data
```matlab
dataSource = 'friedman'
```
Figure 4.

3. Glv simulation data
```matlab
dataSource = 'glv'
```
Figure 2 and 3.

## License
BSD License

## Citation
Hyun-Seob Song, *et al.*,``Minimal Interspecies Interaction Adjustment (MIIA): inference of member-dependent interactions in microbiomes,`` *Bioinformatics*, submitted.

## Contacts
