This folder contains all files that are needed to execute the simulation. The subfolders contain parameter 
files that were used for the simulations in the paper.

Every Simulation can be run by typing the following command in the Windows command line:
Program Parameter_file save_file folder_to_grid

The result is saved as an ASCII file with the x-axis as the first the real part as second and the 
imaginary part as third column.

Example:
LaserIMD LaserIMD//LaserIMD_Doublet_att12.txt LaserIMD//LaserIMD_Doublet_att12_result.txt grid
for a LaserIMD simulation with att12