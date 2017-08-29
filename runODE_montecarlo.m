function output = runODE_montecarlo(A,sysname,q,saveon)
%run ODE for biomolecular system using reduced models including resource
%sharing for many different parameters. Input the weighted, signed
%adjacency matrix, A. Inputs go in rows after adjacency matrix. External
%inputs go in the last row (n+1)st. If no (n+1)st row, there is no external
%input if using random parameters; activation edges > 0; repression edges
%< 0. Edge order is done by cols first (order that A(:) returns).

%defaults
if nargin < 4
    saveon = false;             %save the output and histogram
    if nargin < 3
        q = 50;                 %number of trials to run
        if nargin < 2
            sysname = 'default';%for output titles
            if nargin < 1
                %adjacency matrix
                A = [0 0;
                    0 0;
                    1 0];
            end
        end
    end
end

%settings
addpath utility parameters
ploton = true;      %create a figure
pfun = @(A) params_dist(A);%parameter function handle
tfinal = 1000;      %simulate from 0 to tfinal
umax = 10;          %sweep through inputs from 0 to umax with a trianglewave input
%init
n = min(size(A));   %number of nodes
yfinal = zeros(n,q);%init output vector
[~,prng] = pfun(A);

%run q simulations with structure of A and parameters drawn from the
%parameter distribution in params_dist
for ii = 1:q
    
    %create function handles: edit this function to change form of the dynamics
    [funs,p] = makefuns(A,pfun);
    
    %setup for ODE
    if all(isfield(p,{'RNAP','Ribo','k1','k2','delta1','delta2','DNA','Kp','K2'}))
        xmax = (p.RNAP*p.Ribo*p.k1*p.k2./(p.delta1*p.delta2)).*p.DNA./(p.Kp.*p.K2);
    else
        xmax = 5*ones(n,1);
    end
    x0 = 10.^(log10(xmax(:)).*rand(n,1));
    f_cnst = @(t,x) dynamicsbio_static(t,x,funs,umax);
    opts = odeset('AbsTol',1e-6,'RelTol',1e-4,'Jacobian',...
        @(t,x) dynamicsbio_jac(t,x,funs,umax));
    
    %run ode
    sol = ode15s(f_cnst,[0 tfinal],x0,opts);
    
    %output
    yfinal(:,ii) = sol.y(:,end);
end

%condition output
yfinal(yfinal < 1e-4) = 0;
yfinal = real(yfinal);

output = struct;
output.yfinal = yfinal;
output.prng = prng;


if ploton
    %plot Monte Carlo results
    fig = figure(8); clf;
    %scatterplot of endpoints
    subplot(211)
    semilogy(real(yfinal'),'.')
    ylim([1e-4,inf])
    xlabel('simulation number')
    ylabel('state concentrations [nM]')
    title([strrep(sysname,'_',' '),', n = ',num2str(n)])
    
    %plot histogram
    subplot(212);
    for jj = 1:n
        histogram(log10(yfinal(jj,:)),'BinMethod','Scott',...
            'Normalization','probability','DisplayStyle','stairs')
        hold on;
    end
    hold off;
    ylabel('Frequency')
    xlabel('log_{10}(steady state concentration) [log_{10}(nM)]')
    
    %save settings and output
    if saveon
        savefile = [strrep(sysname,' ','_'),'_MC_n',num2str(n),'_q',num2str(q)];
        cd figures
        savefig(fig,[savefile,'_MC.fig']);
        cd ..
        cd output
        save([savefile,'.mat'],'yf','prng','A','sysname','q');
        cd ..
        
        close(fig);
    end
end
