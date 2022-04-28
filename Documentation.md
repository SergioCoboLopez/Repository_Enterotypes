---
title: "Stochastic Block Models Reveal a Robust Nested Pattern in Healthy Human Gut Microbiomes"
subtitle: ""
author: Sergio Cobo-L\'opez$^1$, Vinod K. Gupta$^{2,3}$, Jaeyun Sung$^{2,3,*}$, 
Roger Guimer\`a$^{1,4,*}$ and Marta Sales-Pardo$^{1,*}$
csl: "podzemna-voda.csl"
header-includes:
  - \usepackage{color}
  - \usepackage{caption}
  - \usepackage{anysize}
---


This repository is written in c++ and python 2.7.
The code MMSBM.cpp finds latent groups of bacteria and hosts using an Expectation-Maximization algorithm. This code is written in C++.
To compile it, you have to have the library armadillo installed and execute the following instruction:

    g++ -std=c++0x -o MMSBM MMSBM.cpp -larmadillo

Once compiled, to execute the code you have to type

    ./MMSBM Seed K L FoldNumber

Where the first argument is the seed, the second and third correspond to the number of groups of hosts and microbes, and the last argument represents the fold in a cross validation experiment.