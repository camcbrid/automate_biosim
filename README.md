# automate_biosim
Automatically simulate ODE biomolecular models for transcription factor networks in MATLAB

Quick file descriptions:
1) runODE_x.m: run ODE solver on dynamicsbio_x.m
2) dynamicsbio_x.m: contains dynamics
3) makefuns.m: create function handles to be used by dynamicsbio_x.m
4) params_y.m: create relevant parameter struct to be used by makefuns.m

Quick start:
To simulate a system, edit the adjacency matrix, A, in runODE_x.m and either select a random parameter function (params_rand.m or params_dist.m which will automatically give the required number of parameters) or input the desired parameters in params.m. 
The adjacency matrix contains a nonzero number in the (i,j) element if there is an edge from node i to node j. The adjacency matrix may be weighted and uses positive values for activation (e.g. +1) and negative values (e.g. -1) for repression.
Parameters for edges in params.m should be specified in the order that is returned by A(:) (i.e. columns of the adjacency matrix first)
External inputs may be specified by adding additional rows to the adjacency matrix

Running methods:
1) x = 'sweep': sweeps input, u, from 0 to umax (specified in runODE_sweep.m)
2) x = 'static': runs system with constant or no external input and outputs equilibrium point. Initial condition is choosen randomly
3) x = 'montecarlo': runs system q times using randomly choosen parameters each run with constant or no external input and outputs the distribution of the equilibrium points. Initial conditions are choosen randomly

Preconfigured systems:
1) y = 'ActCasc': use with runODE_sweep.m. Shows bifurcation in a two node activation cascade due to protease sharing. A = [0 1; 0 0; 1 0]
2) y = 'RepCasc': use with runODE_sweep.m. Shows bifurcation in a two node repression cascade due to ribosome sharing. A = [0 -1; 0 0; -1 0]
3) y = 'ARosc': use with runODE_static.m (with no external input). Shows oscillations in the Hasty oscillator topology. A = [1 1; -1 -1]
