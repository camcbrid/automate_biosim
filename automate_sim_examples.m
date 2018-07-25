
addpath utility parameters


%examples
%activator/repressor oscillator
A = [-1 -1;
    1 1];
pfunARosc = @(A) paramsARosc(A);
outputARosc = runODE_static(A,pfunARosc,0);

%repression cascade with resource sharing
A = [0 -1;
    0 0;
    -1 0];
pfunRepCasc = @(A) paramsRepCasc(A);
urng = [86,94];
outputRepCasc = runODE_sweep(A, pfunRepCasc, urng);

%
A = round(1.5*(rand(5) - .5));
pfundist = params_dist(A);
runODE_static(A,pfundist,0,true)