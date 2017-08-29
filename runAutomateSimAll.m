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
% A = [0 1 -1;
%     0 0 1;
%     -1 0 0];
% q = 100;
% sysname = 'IFFL with neg feedback';

A = [0 0 1 1 1;
    0 0 1 1 1;
    0 0 0 0 0;
    0 0 0 0 0;
    0 0 0 0 0];
n = min(size(A));
q = 100;
sysname = 'Dense overlap pos forward';
pfun = @params_dist;
u = 10;
urng = [0,u];

disp('fmincon')
outfmincon = runfminconbio(A,pfun,u);
disp('fsolve')
outfsolve = runfsolvebio(A,pfun,u);
disp('ODE_static')
outstatic = runODE_static(A,pfun,u);
disp('ODE_sweep')
outsweep = runODE_sweep(A,pfun,urng);
disp('ODE_mulambda')
outmulambda = runODE_mulambda(A,pfun,u);
disp('Monte Carlo')
outMC = runODE_montecarlo(A,sysname,q,false);

%output
dataout = struct;
dataout.fmincon = outfmincon;
dataout.static = outstatic;
dataout.sweep = outsweep;
dataout.fsolve = outfsolve;
dataout.mulambda = outmulambda;
dataout.MC = outMC;

cd output
filename = [sysname,'_n',num2str(n),'_allsimout'];
save([filename,'.mat'],'dataout','A','u','urng','pfun','q','sysname')
cd figs
for ii = 1:8
    if ishandle(ii)
        fig = figure(ii);
        savefig(fig,[filename,num2str(ii),'.fig'])
    end
end
cd ..
cd ..
