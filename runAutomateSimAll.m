%function dataout = runAutomateSimAll(A,u,pfun,urng,sysname,q)
%run each of the possible simulation functions.
%fsolve: find one equilibrium point
%fmincon: find multiple equilibria using fmincon and MultiSearch
%ODE_static: simulate system with constant external input
%ODE_sweep: simulate system sweeping external input with a range
%ODE_mulambda: simulate system varying resource sharing terms
%ODE_montecarlo: run Monte Carlo simulation for parameter distribution to
%find eq pt distribution

%defaults
A = [0 1 -1;
    0 0 1;
    -1 0 0];
q = 1000;
sysname = 'IFFL with neg feedback';
pfun = @params_dist;
u = 10;
urng = [0,u];

disp('fmincon')
outfmincon = runfminconbio(A,pfun,u);
return
disp('fsolve')
outfsolve = runfsolvebio(A,pfun,u);
disp('ODE_static')
outstatic = runODE_static(A,pfun,u);
disp('ODE_sweep')
outsweep = runODE_sweep(A,pfun,urng);
disp('ODE_mulambda')
outmulambda = runODE_mulambda(A,pfun,u);
disp('Monte Carlo')
outMC = runODE_montecarlo(A,sysname,q,true);

%output
dataout = struct;
dataout.fmincon = outfmincon;
dataout.static = outstatic;
dataout.sweep = outsweep;
dataout.fsolve = outfsolve;
dataout.mulambda = outmulambda;
dataout.MC = outMC;
