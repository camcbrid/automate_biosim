%run monte carlo simulations for different graph structures

q = 100;    %number of simulations to run
%A is the adjacency matrix

runODE_montecarlo([1],'positive feedback',q);

n = 10; %number of nodes in the random Erdos Renyi graph
A1 = [rand(n) < 0.8;1,zeros(1,n-1)];
runODE_montecarlo(A1,'ER-0.8',q);

A2 = [rand(n) < 0.5;1,zeros(1,n-1)];
runODE_montecarlo(A2,'ER-0.5',q);

A3 = [rand(n) < 0.3;1,zeros(1,n-1)];
runODE_montecarlo(A3,'ER-0.3',q);

A4 = [rand(n) < 0.1;1,zeros(1,n-1)];
runODE_montecarlo(A4,'ER-0.1',q);

A5 = [0 1 0 0 0;
    0 0 1 0 0;
    0 0 0 1 0;
    0 0 0 0 1;
    1 0 0 0 0;
    1 0 0 0 0];
runODE_montecarlo(A5,'positive feedback',q);

A6 = [0 1 1 1 1;
    1 0 1 1 1;
    1 1 0 1 1;
    1 1 1 0 1;
    1 1 1 1 0;
    1 0 0 0 0];
runODE_montecarlo(A6,'complete positive',q);

A7 = [0 1 0 0 0;
    0 0 1 0 0;
    0 0 0 1 0;
    0 0 0 0 1;
    -1 0 0 0 0;
    1 0 0 0 0];
runODE_montecarlo(A7,'negative feedback',q);

A8 = [0 -1 -1 -1 -1;
    -1 0 -1 -1 -1;
    -1 -1 0 -1 -1;
    -1 -1 -1 0 -1;
    -1 -1 -1 -1 0;
    1 0 0 0 0];
runODE_montecarlo(A8,'complete negative odd',q);

A9 = [0 -1 -1 -1 -1 -1;
    -1 0 -1 -1 -1 -1;
    -1 -1 0 -1 -1 -1;
    -1 -1 -1 0 -1 -1;
    -1 -1 -1 -1 0 -1;
    -1 -1 -1 -1 -1 0
    1 0 0 0 0 0];
runODE_montecarlo(A9,'complete negative even',q);

