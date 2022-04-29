This repository is written in c++ and python 2.7.
The code MMSBM.cpp finds latent groups of bacteria and hosts using an Expectation-Maximization algorithm. This code is written in C++.
To compile it, you have to have the library armadillo installed and execute the following instruction:

    g++ -std=c++0x -o MMSBM MMSBM.cpp -larmadillo

Once compiled, to execute the code you have to type

    ./MMSBM Seed K L FoldNumber

Where the first argument is the seed, the second and third correspond to the number of groups of hosts and microbes, and the last argument represents the fold in a cross validation experiment.

You need to specify input and output paths.
The code takes a training dataset as an input. This dataset is specified by the patient/host number and the fold number. Since we performed a 5-fold cross validation for each patient, there are 5 available folds for each patient.
The code generates three outputs:
1. The loglikelihood file
2. The parameters file: values of etas and theta vectors for each host and microbes, and p matrices.
3. The scores file: how does the algorithm perform at predicting non observed abundances.