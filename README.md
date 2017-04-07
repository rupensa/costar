# CoStar

CO-clustering for STAR-structured data (CoStar) is a multi-view co-clustering algorithm.
Given a set of different views of the same dataset CoStar provides, at the same time, a clustering on documents
and a clustering on features of each view.

If using the code, please cite the following paper:

[1] Ienco, Dino, et al., 2013. Parameter-less co-clustering for star-structured heterogeneous data. Data Mining and Knowledge Discovery 26.2: 217-254

## Run the example

* Compile the sources 

	make all 

* Enter the example data folder 

	cd execution-example

* Run the program with the example data contained in the folder

	../multiViewCoClu -r 6 -n 100 -i config.txt -k 2 -o output

## Parameters of the executable

The executable produced after sources compilation takes in input the following parameters: 

* `-r <int>`, the number of rows (documents) contained in the dataset
* `-n <int>`, the number of costar iterations to perform. Consider that each iteration corresponds to a single element move. 
* `-i <filepath>`, the path (or simply the name) of the configuration file of the current dataset
* `-o <string>`, prefix for the output files
* `-k <int>`, the minimum number of cluster to give as output

### Dataset files

The dataset should be composed by: 

* a txt file for each considered view, these files contains a row for each document in the dataset. 
Each row is a space separated sequence of values
* a configuration file, that contains a row for each view in the dataset. 
Each row contains two space separated values: <n_features, view_file>, where
  * `n_features` is the number of features of the considered view, and 
  * `view_file` is the name of the file containing values for the view.

### Output files

The program produces the following files, where 'prefix' is the prefix specified from the command line.

* prefix.row - contains the row clusters assignment
* prefix_N.col - N files with N from 0 to n_views. Each file contains the column clusters assignment for the i-th view.
* prefix.time - contains the execution time
* prefix.iter - contains for each iteration the tau_r and tau_c values

## Contributions

Original implementation by Dino Ienco, refactoring by [Valentina Rho](https://github.com/valentinarho).

Random Number generator library of type 'Mersenne Twister', by Agner Fog (http://www.agner.org/random/) described in:

M. Matsumoto & T. Nishimura. Mersenne Twister: A 623-Dimensionally Equidistributed Uniform Pseudo-Random Number Generator.
ACM Transactions on Modeling and Computer Simulation, Vol. 8, no. 1, 1998, pp. 3-30.

## Licence

GNU General Public License v3.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

For details of the GNU General Public License refer to <http://www.gnu.org/licenses/>. 
