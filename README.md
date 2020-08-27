# Tesi
 
It is included the Makefile to compile the files.
creation contains a c program that creates a random binary matrix with the increased probability of printing 1 than 0 (80% of 1 and 20% of 0). It can be changed to generate other input structures.

execute creation with "./creation > filename.in"

execute CCL application with "./omp-CCL < input/filename.in > filename.out"
The application prints "Correct" if the result effectively sets correct separations between components. If errors occure, the last wrong value and its position it is displayed.
