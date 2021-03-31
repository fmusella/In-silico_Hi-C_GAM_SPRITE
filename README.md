# *In-silico* implementation of Hi-C, GAM and SPRITE methods
## Authors
Luca Fiorillo and Francesco Musella
## Reference
Fiorillo et al., *“Comparison of the Hi-C, GAM and SPRITE methods by using polymer models of chromatin”* (2021).
## Details
We provide the algorithms to simulate simplified versions of Hi-C, GAM and SPRITE experiments over model polymer structures. The algorithms require as input the spatial coordinates of polymer 3D structures. To account for the two alleles present in real single cells, we consider pairs of polymer structures, which we refer to as *in-silico* cells. By iterating the algorithms over the *in-silico* cells we simulate experiments carried over a population of cells. Here, we provide the SBS model 3D conformations of the Sox9 locus (mESC chr11:109-115Mb, mm9). The output of the Hi-C and SPRITE algorithms is an *in-silico* contact matrix, showing the detected contacts between pairs of polymer segments across *in-silico* cells. For the GAM algorithm, the output is an *in-silico* co-segregation matrix, reporting the frequencies whereby pairs of polymer segments are detected in random slices extracted from the *in-silico* cells. 

All codes are written in C language. The source codes and their headers are respectively in the *src* and *headers* folders. The executable files are in the *main* folder. In the setting.c file and its header, parameters relative to the simulations are set. In the *data* folder, the *structures* subfolder contains the model 3D polymer structures of the Sox9 locus; the *pairs* subfolder contains pairs of structures arranged in a spherical volume so that a random slice can be extracted from them in the GAM algorithm. 
## How to run
To run the computational implementations of Hi-C, SPRITE and GAM over SBS polymer structures the following steps are needed:
1. In the *main* directory, compile the C program main.c using the customized compiler compile.sh with the following command:
```
./compile.sh main.c
```
2. That generates the executable main.out. It must be launched on the command line according to the following syntax:
```
./main.out   which_method   Ncells   efficiency
```
  -	which_method defines which one of the Hi-C, SPRITE and GAM simulation algorithms will be run. It must be 1 to run the Hi-C algorithm, 2 to run the SPRITE algorithm and 3 to run the GAM algorithm. 
  -	Ncells indicates the number of in-silico cells to employ.
  -	efficiency indicates the efficiency whereby the computational version of Hi-C, SPRITE or GAM is performed.

For example, to run the in-silico version of Hi-C over 100 *in-silico* cells at 0.05 efficiency, the following command must be typed:
```
./main.out 1 100 0.05
```
The codes were tested on macOS High Sierra 10.13.6 and on Windows 10 (using MobaXterm_Personal 20.1). 
