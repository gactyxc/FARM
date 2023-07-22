# FARM
the source code of FARM method for automatic rice mapping.
# Code introduction
## step 1: get the SAR features 
It should be noticed that if change the region, the s1 features should be firstly exported to your Asset using the export code.
## step 2: get the land cover objects based on SNIC
## step 3 (option step): get the T1 and T2 value for sample selection
A supply code for getting the mor accurate values of T1 and T2 in FARM. 
The default value of T1 was -20; the default value of T2 was -17 in central and southern China and -20 fro northern China.
## step 4: get the rice and non-rice objects 
## step 5: classification using the multi-RF
