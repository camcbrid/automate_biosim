function output = runLHS_static(A, numsamples, u)
%OUTPUT = runLHS_static(A, NUMSAMPLES, U)
%run multiple simulations using Latin Hypercube sampling to sample
%parameter space. Input adjacency matrix, number of samples, and external
%input concentration. Uses runODE_static. A is a weighted adjacency matrix,
%NUMSAMPLES is the number of samples to run for the Monte Carlo simulation,
%and U is the input to the circuit (on the (n+1)st row of the adjacency
%matrix.

%defauts
if nargin < 3
    u = 1;
    if nargin < 2
        numsamples = 250;
        if nargin < 1
            %adjacency matrix
            A = [0 0 1 1 1;
                0 0 1 1 1;
                0 0 0 0 0;
                0 0 0 0 0;
                0 0 0 0 0];
        end
    end
end

%parameter range
p = params_LHS(A,numsamples);

%settings
addpath parameters utility
ploton = true;          %if want to display a figure
tfinal = 50;            %simulate from 0 to tfinal
n = min(size(A));       %number of nodes
%init
yf = zeros(n,numsamples);
yout = cell(0);

%loop through each Latin Hypercube sample
for ii = 1:numsamples
    
    tic
    %create function handles: edit this function to change form of the dynamics
    funs = makefuns(A,p{ii});
    
    %setup for ODE
    if all(isfield(p{ii},{'RNAP','Ribo','k1','k2','delta1','delta2','DNA','Kp','K2'}))
        xmax = (p{ii}.RNAP*p{ii}.Ribo*p{ii}.k1*p{ii}.k2./...
            (p{ii}.delta1*p{ii}.delta2)).*p{ii}.DNA./(p{ii}.Kp.*p{ii}.K2);
    else; xmax = 5*ones(n,1);
    end
    x0 = 10.^(log10(xmax(:)).*rand(n,1));
    f_cnst = @(t,x) dynamicsbio_static(t,x,funs,u,1,1);
    opts = odeset('AbsTol',1e-6,'RelTol',1e-4,'Stats','off','Vectorized','on',...
        'Jacobian',@(t,x) dynamicsbio_jac(t,x,funs,u));
    
    %run ode
    sol = ode15s(f_cnst,[0 tfinal],x0,opts);
    %outputs
    tout = linspace(0,tfinal,700)';
    yout{ii} = deval(sol,tout);
    yf(:,ii) = yout{ii}(:,end);
    toc
end

%output data
output = struct;
output.t = tout;
output.y = yout;
output.yf = yf;
output.p = p;

%plot
if ploton
    %plot output
    figure(9); clf;
    %final values
    subplot(211)
    semilogy(yf')
    title('Latin Hypercube Distribution')
    ylabel('tfinal protein concentration [nM]')
    xlabel('sample')
    %histogram of final values
    subplot(212)
    histogram(log10(yf'),'BinMethod','Scott',...
        'Normalization','probability','DisplayStyle','stairs');
    ylabel('frequency')
    xlabel('log_{10}(steady state concentration) [log_{10}(nM)]')
end