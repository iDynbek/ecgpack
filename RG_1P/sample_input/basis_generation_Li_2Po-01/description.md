## Description

Sample file `inout.txt` in this directory provides input for generating ECG basis for the lowest $^2 P^o$ state of Li atom (Li-7 isotope). It starts from 0 basis functions and stops when the basis reaches 100 basis functions. When the basis size equals 5, 10, 20, 50, and 100 functions the program will save the input in separate files named  
`inout_Li_2Po-01-00005.txt`  
`inout_Li_2Po-01-00010.txt`  
`inout_Li_2Po-01-00020.txt`  
`inout_Li_2Po-01-00050.txt`  
`inout_Li_2Po-01-00100.txt`  
 The calculation with this input file should take roughly 1-2 minutes on a single CPU core when double precision is used (gfortran compiler / AMD Ryzen 7 7800X3D).
